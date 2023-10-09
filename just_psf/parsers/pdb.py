import numpy

from just_psf import logger
from just_psf.geometry import PDBGeometry
from just_psf.parsers import ParseError
from just_psf.parsers.line import TokenType, LineParser


l_logger = logger.getChild(__name__)


class PDBParseError(ParseError):
    pass


class PDBParser(LineParser):
    """Parse a Protein Data Bank (PDB) file.

    Format is defined at https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html.

    Note: actually only parse `ATOM`/`HETATM`/`END`!
    """

    def pdb(self):
        seg_names = []
        resi_ids = []
        resi_names = []
        atom_types = []
        positions = []
        symbols = []

        while self.current_token.type == TokenType.LINE:
            kw = self.current_token.value[:6].strip()
            if kw == 'END':
                break

            elif kw in ['ATOM', 'HETATM']:
                li = self.current_token.value
                self.next()

                if len(li) < 78:
                    raise PDBParseError(self.current_token, 'len(line) < 78')

                serial, atype, resi_name, seg_name, res_id, x, y, z, symbol = (
                    li[6:11].strip(),  # serial
                    li[12:16].strip(),  # name
                    # li[16].strip(),  # altloc
                    li[17:21].strip(),  # resname
                    li[21].strip(),  # chainid
                    li[22:26].strip(),  # resseq
                    # li[26],  # icode
                    li[30:38].strip(),  # x
                    li[38:46].strip(),  # y
                    li[46:54].strip(),  # z
                    # li[54:60].strip(),  # occ
                    # li[60:66].strip(),  # tempfac
                    # li[66:76].strip(),  # segID
                    li[76:78].strip(),  # element
                )

                if serial == '':
                    raise PDBParseError(self.current_token, 'no serial')
                if symbol == '':
                    raise PDBParseError(self.current_token, 'no symbol')
                if '' in [x, y, z]:
                    raise PDBParseError(self.current_token, 'empty coordinates')

                seg_names.append(seg_name)
                resi_ids.append(int(res_id))
                resi_names.append(resi_name)
                atom_types.append(atype)
                symbols.append(symbol)
                positions.append(tuple(float(n) for n in [x, y, z]))

            else:
                self.next()  # skip that line

            self.next_non_empty()

        if not self.current_token.value == 'END':
            raise PDBParseError(self.current_token, 'expected `END`')

        self.eat(TokenType.LINE)
        self.next_non_empty()
        self.eat(TokenType.EOF)

        return PDBGeometry(
            symbols=symbols,
            positions=numpy.array(positions).reshape(-1, 3),
            seg_names=seg_names,
            resi_ids=resi_ids,
            resi_names=resi_names,
            atom_types=atom_types
        )
