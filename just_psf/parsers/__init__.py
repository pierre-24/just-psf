
# parsers error
class ParseError(Exception):
    def __init__(self, token, msg: str):
        super().__init__('on line {}: {}'.format(token.line, msg))
        self.token = token
