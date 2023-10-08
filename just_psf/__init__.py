"""
Just get me a topology, for god’s sake!
"""

import logging
import os

__version__ = '0.1.0'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


# logging
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
logger.setLevel(os.environ.get('LOGLEVEL', 'WARNING').upper())


# error
class ParseError(Exception):
    def __init__(self, token, msg: str):
        super().__init__('on line {}: {}'.format(token.line, msg))
        self.token = token
