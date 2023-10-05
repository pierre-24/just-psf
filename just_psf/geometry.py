import numpy

from typing import TextIO
from just_psf import logger


l_logger = logger.getChild(__name__)


class Geometry:
    def __init__(self, symbols: list, positions: numpy.ndarray):
        assert positions.shape == (len(symbols), 3)

        self.symbols = symbols
        self.positions = positions

    def __len__(self) -> int:
        return len(self.symbols)

    def copy(self) -> 'Geometry':
        """Copy itself. Involves a copy of positions and symbols.
        """

        return Geometry(
            self.symbols.copy(),
            self.positions.copy()
        )

    @classmethod
    def from_xyz(cls, f: TextIO) -> 'Geometry':
        """Read geometry from a XYZ file
        """

        l_logger.debug('Reading geometry...')

        symbols = []
        positions = []

        n = int(f.readline())
        f.readline()

        for i in range(n):
            data = f.readline().split()
            symbols.append(data[0])
            positions.append([float(x) for x in data[1:]])

        l_logger.debug('... Got {} atom(s)'.format(n))

        return cls(symbols, numpy.array(positions))

    def to_xyz(self, title: str = '') -> str:
        """Get XYZ representation of this geometry"""

        r = '{}\n{}'.format(len(self), title)
        for i in range(len(self)):
            r += '\n{:2} {: .7f} {: .7f} {: .7f}'.format(self.symbols[i], *self.positions[i])

        return r
