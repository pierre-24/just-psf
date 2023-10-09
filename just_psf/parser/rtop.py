import io
from typing import Iterator, TextIO, Optional, Callable, List, Tuple, Set, Union
from enum import Enum, unique

import numpy

from just_psf import logger
from just_psf.parser import ParseError

from just_psf.residue_topology import ResidueTopology, Topologies


l_logger = logger.getChild(__name__)


@unique
class TokenType(Enum):
    WORD = 'WRD'
    TITLEL = 'TIL'
    NL = 'NWL'
    EOF = '\0'


COMMENT = '!'
TITLEL = '*'
SPACES = [' ', '\t']
NLS = ['\n', '\r']


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


class RTopParseError(ParseError):
    pass


class RTopParser:
    """Parse a RTop/RTF (Residue Topology(ies) File) file from CHARMM.
    Also referred to as "toppar" (in CHARMM files).

    Do not parse `GROU`, `IC`, `IMPR`, `CMAP`, `DONO`, `ACCE`, `PATC`, but just skip them.

    Sources:
    + http://www.ks.uiuc.edu/Training/Tutorials/science/topology/topology-html/node4.html
    + https://www.charmm-gui.org/?doc=lecture&module=molecules_and_topology&lesson=2
    """

    def __init__(self, inp: Union[TextIO, str]):
        if isinstance(inp, io.TextIOBase):
            self.input = inp.read()  # TODO: is it possible to stream the result like PSF?
        else:
            self.input = inp

        self.current_token: Optional[Token] = None
        self.current_line = 1
        self.position = 0

        self.next()

    def tokenize(self) -> Iterator[Token]:
        while self.position < len(self.input):
            start = self.position
            line_start = self.current_line

            if self.input[start] in SPACES:  # skip spaces
                self.position = self._skip_if(predicate=lambda x: x in SPACES)
            elif self.input[start] in NLS:
                self.position += 1
                self.current_line += 1
                yield Token(TokenType.NL, self.input[start:self.position], line_start)
            else:
                if self.input[start] == '*' and (start == 0 or self.input[start - 1] in NLS):
                    self.position = self._skip_if(predicate=lambda x: x not in NLS)
                    yield Token(TokenType.TITLEL, self.input[start:self.position], line_start)
                elif self.input[start] == COMMENT:  # skip comment
                    self.position = self._skip_if(predicate=lambda x: x not in NLS)
                else:
                    self.position = self._skip_if(predicate=lambda x: x not in SPACES and x not in NLS)
                    yield Token(TokenType.WORD, self.input[start:self.position], line_start)

        yield Token(TokenType.EOF, '\0', self.position)

    def _skip_if(self, predicate: Callable) -> int:
        """Go to the next position while `predicate` is true"""
        end = self.position + 1
        while end < len(self.input) and predicate(self.input[end]):
            end += 1

        return end

    def next(self):
        self.current_token = next(self.tokenize())

    def next_non_empty(self):
        """Got to the next non-NL token
        """
        while self.current_token.type == TokenType.NL:
            self.next()

    def expect(self, typ: TokenType):
        if self.current_token.type != typ:
            raise RTopParseError(self.current_token, 'expected {}, got {}'.format(typ, self.current_token))

    def eat(self, typ: TokenType):
        if self.current_token.type == typ:
            self.next()
        else:
            raise RTopParseError(self.current_token, 'expected {}, got {}'.format(typ, self.current_token))

    def word(self) -> str:
        """Parse word
        """

        self.expect(TokenType.WORD)
        word = self.current_token.value
        self.next()
        return word

    def integer(self) -> int:
        """Parse integer
        """
        self.expect(TokenType.WORD)

        try:
            number = int(self.current_token.value)
        except ValueError:
            raise RTopParseError(self.current_token, 'expected integer, got {}'.format(self.current_token.value))

        self.next()
        return number

    def number(self) -> float:
        """Parse float
        """
        self.expect(TokenType.WORD)

        try:
            number = float(self.current_token.value)
        except ValueError:
            raise RTopParseError(self.current_token, 'expected number, got {}'.format(self.current_token.value))

        self.next()
        return number

    def topologies(self) -> Topologies:
        """
        TOPOLOGY := HEADER DECLS RESIDUE* END
        HEADER := TITLE VERSION (MASS | DECL | DEFA | AUTO)*
        MASS := 'MASS' INT STRING FLOAT WORD? NL
        DECL := 'DECL' WORD NL
        DEFA := 'DEFA' (STRING STRING)* NL
        AUTO := 'AUTO' WORD* NL
        """

        l_logger.debug('Parsing topologies')

        top_masses = {}
        top_autogenerate = set()
        top_decls = []
        top_defaults = set()

        # header
        self.title()
        major, minor = self.version()
        l_logger.debug('Rtop v{}.{}'.format(major, minor))

        self.next_non_empty()

        # a few declarations before the residues
        while self.current_token.type == TokenType.WORD and self.current_token.value[:4] not in ['RESI', 'END']:
            keyword = self.current_token.value[:4]
            self.next()
            if keyword == 'MASS':
                self.integer()  # skip id

                atyp = self.word()
                if atyp in top_masses:
                    raise RTopParseError(self.current_token, 'a MASS is already defined for `{}`'.format(atyp))

                mass = self.number()

                if self.current_token.type == TokenType.WORD:  # skip mmff
                    self.next()

                top_masses[atyp] = mass
            elif keyword == 'AUTO':
                autogen = []
                while self.current_token.type == TokenType.WORD:
                    autogen.append(self.word())
                top_autogenerate.add(tuple(autogen))
            elif keyword == 'DECL':
                w = self.word()
                if w in top_decls:
                    raise RTopParseError(self.current_token, '`{}` is already DECLared'.format(w))
                top_decls.append(w)
            elif keyword == 'DEFA':
                while self.current_token.type == TokenType.WORD:
                    top_defaults.add((self.word(), self.word()))
            else:
                raise RTopParseError(self.current_token, 'unknown keyword `{}`'.format(keyword))

            if self.current_token.type != TokenType.EOF:
                self.eat(TokenType.NL)
                self.next_non_empty()

        # residues
        allowed_types = set(top_masses.keys())
        residues = []
        while self.current_token.type == TokenType.WORD and self.current_token.value[:4] == 'RESI':
            residues.append(self.residue(allowed_types, top_decls))

        # normally, there is nothing more:
        end_keyword = self.word()
        if end_keyword != 'END':
            raise RTopParseError(self.current_token, 'expected `END`, got `{}`'.format(end_keyword))

        self.next_non_empty()
        self.eat(TokenType.EOF)

        l_logger.debug('Done')

        return Topologies(
            masses=top_masses,
            autogenerate=top_autogenerate,
            defaults=top_defaults,
            declarations=top_decls,
            residues=residues
        )

    def title(self) -> List[str]:
        """
        TITLE := ((TITLEL NL)* '*' NL)?
        """

        title_lines = []

        if self.current_token.type == TokenType.TITLEL:
            while self.current_token.type == TokenType.TITLEL and self.current_token.value != '*':
                title_lines.append(self.current_token.value[1:])
                self.next()
                self.eat(TokenType.NL)

            if not self.current_token.type == TokenType.TITLEL:
                raise RTopParseError(self.current_token, 'expected `*` to end the title')

            self.next()
            self.eat(TokenType.NL)

        return title_lines

    def version(self) -> Tuple[int, int]:
        """
        VERSION := INT INT NL
        """

        major_version = self.integer()

        if major_version < 19:
            raise RTopParseError(self.current_token, 'cannot parse RTop file version < 19')

        minor_version = self.integer()
        self.eat(TokenType.NL)

        return major_version, minor_version

    def residue(self, allowed_types: Set[str], declarations: List[str]) -> ResidueTopology:
        """
        RESIDUE := 'RESI' WORD NUMBER NL DEF*
        DEF := ATOM | GROUP | BOND | DOUBLE | IMPR | CMAP | DONOR | ACCEPTOR | IC
        """

        atom_names = []
        atom_types = []
        atom_charges = []
        bonds = []

        name_to_index = {}
        for i, n in enumerate(declarations):
            name_to_index[n] = -1 - i

        if self.current_token.type != TokenType.WORD or self.current_token.value[:4] != 'RESI':
            raise RTopParseError(self.current_token, 'expected `RESI`')

        self.next()
        name = self.word()
        charge = self.number()
        self.eat(TokenType.NL)

        l_logger.debug('Parsing residue `{}`'.format(name))

        while self.current_token.type == TokenType.WORD and self.current_token.value[:4] not in ['RESI', 'END']:
            keyword = self.current_token.value[:4]
            self.next()

            if keyword == 'ATOM':
                aname = self.word()
                if aname in name_to_index:
                    raise Exception('an atom with name `{}` already exists'.format(aname))
                name_to_index[aname] = len(atom_names)
                atom_names.append(aname)

                atype = self.word()
                if atype not in allowed_types:
                    raise RTopParseError(self.current_token, 'unknown atom type `{}`'.format(atype))
                atom_types.append(atype)

                atom_charges.append(self.number())

            elif keyword in ['BOND', 'DOUB']:
                while self.current_token.type == TokenType.WORD:
                    a1, a2 = self.word(), self.word()
                    try:
                        bonds.append((name_to_index[a1], name_to_index[a2]))
                    except KeyError as e:
                        raise Exception('unknown atom name `{}` in bond'.format(e))

            elif keyword in ['GROU', 'IC', 'IMPR', 'CMAP', 'DONO', 'ACCE', 'PATC']:  # just skip
                while self.current_token.type == TokenType.WORD:
                    self.next()
            else:
                raise RTopParseError(self.current_token, 'unknown keyword `{}` in RESI'.format(keyword))

            if self.current_token.type != TokenType.EOF:
                self.eat(TokenType.NL)
                self.next_non_empty()

        return ResidueTopology(
            resi_name=name,
            resi_charge=charge,
            atom_names=atom_names,
            atom_types=atom_types,
            atom_charges=atom_charges,
            bonds=numpy.array(bonds).reshape((-1, 2))
        )
