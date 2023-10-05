from typing import TextIO, Iterable, Optional

from enum import Enum, unique
import numpy
from numpy.typing import NDArray

from just_psf import logger
from just_psf.topology import Topology


l_logger = logger.getChild(__name__)


class PSFParseError(Exception):
    pass


@unique
class TokenType(Enum):
    LINE = 'LIN'
    EMPTY = 'EMP'
    EOF = 'EOF'


class Token:
    def __init__(self, type: TokenType, value: str, line: int = -1):
        self.type = type
        self.value = value
        self.line = line

    def __repr__(self):
        return '<Token({},{}{})>'.format(
            self.type,
            repr(self.value),
            '' if self.line < 0 else ', {}'.format(self.line)
        )


class PSFParser:
    """Parse a PSF file.
    Tries to follow as closely as possible the actual PSF format (e.g., fixed length for fields).
    Handle the `NAMD` (but only for the atom section) `EXT` and `XPLOR` flags.
    Does not parse `CHEQ` data (but should not fail to parse the rest!).

    Sources:
    - https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html (quick explanation)
    - https://www.charmm-gui.org/?doc=lecture&module=pdb&lesson=6 (more explanation, thanks people!)
    - https://docs.mdanalysis.org/2.4.1/_modules/MDAnalysis/topology/PSFParser.html#PSFParser (example of read)
    - https://parmed.github.io/ParmEd/html/_modules/parmed/charmm/psf.html#CharmmPsfFile (another example)
    """

    PARSED_SECTIONS = {
        # 'title': (indices_per, elements_per)
        'NATOM': None,
        'NBOND': (2, 4),
        'NTHETA': (3, 3),
        'NPHI': (4, 2),
        'NIMPHI': (4, 2),
        'NDON': (2, 4),
        'NACC': (2, 4)
    }

    def __init__(self, f: TextIO):
        self.source = f
        self.current_token: Optional[Token] = None

        self.next()

    def tokenize(self) -> Iterable[str]:

        n = 1

        while True:
            line = self.source.readline()
            if line == '':
                break
            else:
                n += 1
                if line == '\n':
                    yield Token(TokenType.EMPTY, '', n)
                else:
                    yield Token(TokenType.LINE, line[:-1], n)

        yield Token(TokenType.EOF, '\0')

    def next(self):
        self.current_token = next(self.tokenize())

    def next_if_empty_or_raises(self):
        if self.current_token.type not in [TokenType.EMPTY, TokenType.EOF]:
            raise PSFParseError('on line {}: expected empty line'.format(self.current_token.line))

        self.next()

    def topology(self) -> Topology:
        l_logger.debug('Parsing topology')

        sections = {}
        atoms = None

        if not self.current_token.value.startswith('PSF'):
            raise PSFParseError('this is not a PSF file')

        flags = self.current_token.value.split()[1:]
        self.next()
        self.next_if_empty_or_raises()

        intsize = 10 if 'EXT' in flags else 8
        l_logger.debug('Will use intsize={}'.format(intsize))

        while self.current_token.type != TokenType.EOF:
            # read section:
            pos = self.current_token.value.find('!')
            if pos < 0:
                raise PSFParseError('on line {}: expected section'.format(self.current_token.line))

            title = self.current_token.value[pos + 1:]
            if ':' in title:
                title = title[:title.index(':')]

            l_logger.debug('Found section `{}`'.format(title))

            if title in self.PARSED_SECTIONS:
                l_logger.debug('Parsing...')
                try:
                    n = int(self.current_token.value[:intsize])
                except ValueError:
                    raise PSFParseError('on line {}: incorrectly formatted section'.format(self.current_token.line))

                self.next()

                if title == 'NATOM':
                    atoms = self.parse_atoms(n, flags)
                else:
                    sections[title] = self.parse_indices(intsize, n, *self.PARSED_SECTIONS[title])

                l_logger.debug('... Got {} elements'.format(n))

                self.next_if_empty_or_raises()  # a section ends with an empty line
            else:
                self.next()
                self.skip_section()
                l_logger.debug('Skipped')

        if atoms is None:
            raise PSFParseError('missing NATOM section')

        return Topology(
            *atoms,
            sections.get('NBOND', None),
            sections.get('NTHETA', None),
            sections.get('NPHI', None),
            sections.get('NIMPHI', None),
            sections.get('NDON', None),
            sections.get('NACC', None)
        )

    def skip_section(self):
        """Skip to the next section
        """

        while self.current_token.type != TokenType.EOF and '!' not in self.current_token.value:
            self.next()

    def parse_atoms(self, n: int, flags: list):
        """Inspired by
        https://docs.mdanalysis.org/2.4.1/_modules/MDAnalysis/topology/PSFParser.html#PSFParser (without CHEQ!)
        Format checked against `charmm/source/io/psfres.F90` (in c47b1).
        """

        atom_ids = []
        seg_names = []
        resi_ids = []
        resi_names = []
        symbols = []
        atom_types = []
        charges = []
        masses = []
        fixed = []

        atom_parsers = {
            'STANDARD': lambda li:
            # fmt02='(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)' (XPLOR use A4 for atom_type instead of I4)
            (
                li[:8],  # I8 atom_id
                # 1X
                li[9:13].strip() or 'SYS',  # A4 seg_name
                # 1X
                li[14:18],  # A4 resi_id
                # 1X
                li[19:23].strip(),  # A4 resi_name
                # 1X
                li[24:28].strip(),  # A4 symbols
                # 1X
                li[29:33].strip(),  # A4  atom_type
                # 1X
                li[34:48], li[48:62],  # 2G14.6 charge mass
                li[62:70]),  # I8 fixed
            'EXTENDED': lambda li:
            # fmt01='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)'
            (
                li[:10],  # I10 atom_id
                # 1X
                li[11:19].strip() or 'SYS',  # A8 seg_name
                # 1X
                li[20:28],  # A8 resi_id
                # 1X
                li[29:37].strip(),  # A8 resi_name
                # 1X
                li[38:46].strip(),  # A8 symbol
                # 1X
                li[47:51].strip(),  # A4 atom_type
                # 1X
                li[52:66], li[66:79],  # 2G14.6 charge mass
                li[79:87]  # I8 fixed
            ),
            'EXTENDED_XPLOR': lambda li:
            #  fmt02='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)'
            (
                li[:10],  # I10 atom_id
                # 1X
                li[11:19].strip() or 'SYS',  # A8 seg_name
                # 1X
                li[20:28],  # A8 resi_id
                # 1X
                li[29:37].strip(),  # A8 resi_name
                # 1X
                li[38:46].strip(),  # A8 symbol
                # 1X
                li[47:53].strip(),  # A6 atom_type
                # 1X
                li[54:68], li[68:81],  # 2G14.6 charge mass
                li[81:89]  # I8 fixed
            ),
            'NAMD': lambda li: li.split()[:8],
        }

        if 'NAMD' in flags:
            l_logger.debug('will use NAMD parser')
            parser = atom_parsers['NAMD']
        elif 'EXT' in flags:
            if 'XPLOR' in flags:
                l_logger.debug('will use EXTENDED_XPLOR parser')
                parser = atom_parsers['EXTENDED_XPLOR']
            else:
                l_logger.debug('will use EXTENDED parser')
                parser = atom_parsers['EXTENDED']
        else:
            l_logger.debug('will use STANDARD parser')
            parser = atom_parsers['STANDARD']

        i = 0
        while self.current_token.type != TokenType.EOF:
            if self.current_token.type == TokenType.EMPTY:
                break

            if i >= n:
                raise PSFParseError('on line {}: too much atoms, expected {}'.format(self.current_token.line, n))

            values = parser(self.current_token.value)

            atom_ids.append(int(values[0]))
            seg_names.append(values[1])
            resi_ids.append(int(values[2]))
            resi_names.append(values[3])
            symbols.append(values[4])
            atom_types.append(values[5])
            charges.append(float(values[6]))
            masses.append(float(values[7]))
            fixed.append(int(values[8]) == 1)

            i += 1
            self.next()

        if i != n:
            raise PSFParseError('on line {}: not enough atoms, expected {}'.format(self.current_token.line, n))

        return symbols, atom_types, charges, atom_ids, seg_names, resi_ids, resi_names, masses, fixed

    def parse_indices(self, intsize: int, n: int, indices_per: int, elements_per: int) -> Optional[NDArray[int]]:
        """Read out `n` elements, each of witch contains `indices_per` indices.
        There are `elements_per` elements per line, written using `intsize` characters.
        """

        atm_ids = numpy.ndarray((n * indices_per, ), dtype=numpy.int64)
        i = 0
        remaining = n

        if n == 0:
            self.next_if_empty_or_raises()  # 0-length section contains an empty line
            return None

        while self.current_token.type != TokenType.EOF:
            if self.current_token.type == TokenType.EMPTY:
                break

            if remaining == 0:
                raise PSFParseError('on line {}: too much data, got already {}'.format(self.current_token.line, n))

            to_read = min(elements_per * indices_per, remaining * indices_per)

            if len(self.current_token.value) != to_read * intsize:
                raise PSFParseError('on line {}: incorrect number of indices, expected {}'.format(
                    self.current_token.line, to_read))

            try:
                atm_ids[i:i + to_read] = numpy.int64(
                    tuple(self.current_token.value[j * intsize:(j + 1) * intsize] for j in range(to_read)))
            except ValueError:
                raise PSFParseError('on line {}: unable to parse indices'.format(self.current_token.line))

            i += to_read
            remaining -= to_read // indices_per
            self.next()

        if remaining != 0:
            raise PSFParseError('on line {}: not enough data, {} missing'.format(self.current_token.line, remaining))

        return atm_ids.reshape(n, indices_per) - 1  # indices are 0-based in `Topology`
