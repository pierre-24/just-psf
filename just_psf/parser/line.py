from typing import TextIO, Iterator, Optional

from enum import Enum, unique
from just_psf.parser import ParseError


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


class LineParser:
    """Read a file line per line
    """

    def __init__(self, f: TextIO):
        self.source = f
        self.current_token: Optional[Token] = None
        self.current_line = 0

        self.next()

    def tokenize(self) -> Iterator[str]:
        line = self.source.readline()

        while line != '':
            self.current_line += 1
            if line == '\n':
                yield Token(TokenType.EMPTY, '', self.current_line)
            else:
                yield Token(TokenType.LINE, line[:-1], self.current_line)

            line = self.source.readline()

        yield Token(TokenType.EOF, '\0')

    def next(self):
        self.current_token = next(self.tokenize())

    def next_non_empty(self):
        while self.current_token.type == TokenType.EMPTY:
            self.next()

    def expect(self, typ: TokenType):
        if self.current_token.type != typ:
            raise ParseError(self.current_token, 'expected {}, got {}'.format(typ, self.current_token.type))

    def eat(self, typ: TokenType):
        if self.current_token.type == typ:
            self.next()
        else:
            raise ParseError(self.current_token, 'expected {}, got {}'.format(typ, self.current_token.type))
