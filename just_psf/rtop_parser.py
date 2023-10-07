from typing import Iterator, TextIO, Optional, Callable, List, Tuple, Dict
from enum import Enum, unique

from just_psf.topology import Residue


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


class RTopParseError(Exception):
    def __init__(self, token: Token, msg: str):
        super().__init__('on line {}: {}'.format(token.line, msg))
        self.token = token


class RTopParser:
    """Parse a RTop (Residue Topology File) file from CHARMM.

    Sources:
    + http://www.ks.uiuc.edu/Training/Tutorials/science/topology/topology-html/node4.html
    + https://www.charmm-gui.org/?doc=lecture&module=molecules_and_topology&lesson=2
    """

    def __init__(self, f: TextIO):
        self.input = f.read()
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

    def topology(self):
        """
        TOPOLOGY := HEADER DECLS RESIDUE* END
        HEADER := TITLE VERSION
        """

        # header
        _ = self.title()
        self.next_non_empty()
        _ = self.version()
        self.next_non_empty()

        # decls
        _ = self.decls()

        # residues
        residues = []
        while self.current_token.type == TokenType.WORD and self.current_token.value[:4] == 'RESI':
            residues.append(self.residue())

        end_keyword = self.word()
        if end_keyword != 'END':
            raise RTopParseError(self.current_token, 'expected `END`, got `{}`'.format(end_keyword))

        self.eat(TokenType.EOF)

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

    def decls(self) -> Tuple[Dict[str, float], List[str], List[Tuple[str, str]], List[str]]:
        """
        DECLS := (MASS | DECL | DEFA | AUTO)*
        MASS := 'MASS' INT STRING FLOAT WORD? NL
        DECL := 'DECL' WORD NL
        DEFA := 'DEFA' (STRING STRING)* NL
        AUTO := 'AUTO' WORD* NL
        """

        self.expect(TokenType.WORD)

        masses = {}
        autos = set()
        decls = set()
        defas = []

        while self.current_token.type != TokenType.EOF:
            keyword = self.current_token.value[:4]
            if keyword == 'RESI':
                break
            else:
                self.next()
                if keyword == 'MASS':
                    self.integer()  # skip id

                    atyp = self.word()
                    if atyp in masses:
                        raise RTopParseError(self.current_token, 'a MASS is already defined for `{}`'.format(atyp))

                    mass = self.number()

                    if self.current_token.type == TokenType.WORD:  # skip mmff
                        self.next()

                    masses[atyp] = mass
                elif keyword == 'AUTO':
                    while self.current_token.type == TokenType.WORD:
                        autos.add(self.word())
                elif keyword == 'DECL':
                    w = self.word()
                    if w in decls:
                        raise RTopParseError(self.current_token, '`{}` is already DECLared'.format(w))
                    decls.add(w)
                elif keyword == 'DEFA':
                    while self.current_token.type == TokenType.WORD:
                        defas.append((self.word(), self.word()))
                else:
                    raise RTopParseError(self.current_token, 'unknown keyword `{}`'.format(keyword))

                if self.current_token.type != TokenType.EOF:
                    self.eat(TokenType.NL)
                    self.next_non_empty()

        return masses, list(decls), defas, list(autos)

    def residue(self) -> Residue:
        """
        RESIDUE := 'RESI' WORD NUMBER NL DEF*
        DEF := ATOM | GROUP | BOND | DOUBLE | IMPR | CMAP | DONOR | ACCEPTOR | IC
        """

        if not self.current_token.type != TokenType.WORD or self.current_token.value[:4] != 'RESI':
            raise RTopParseError(self.current_token, 'expected `RESI`')

        name = self.word()
        charge = self.number()
        self.eat(TokenType.NL)

        return Residue(name, charge)
