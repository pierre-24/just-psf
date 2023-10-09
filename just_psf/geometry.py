import numpy
from typing import TextIO, List, Optional
from numpy.typing import NDArray

from just_psf import logger


l_logger = logger.getChild(__name__)


class Geometry:
    def __init__(self, symbols: List[str], positions: NDArray[float]):
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


class PDBGeometry(Geometry):
    def __init__(
        self,
        symbols: List[str],
        positions: NDArray[float],
        seg_names: Optional[List[str]] = None,
        resi_ids: Optional[List[int]] = None,
        resi_names: Optional[list[str]] = None,
        atom_types: Optional[List[str]] = None,
    ):
        super().__init__(symbols, positions)

        assert seg_names is None or len(seg_names) == len(symbols)
        assert resi_ids is None or len(resi_ids) == len(symbols)
        assert resi_names is None or len(resi_names) == len(symbols)
        assert atom_types is None or len(atom_types) == len(symbols)

        self.seg_names = seg_names
        self.resi_ids = resi_ids
        self.resi_names = resi_names
        self.atom_types = atom_types

    @classmethod
    def from_pdb(cls, f: TextIO) -> 'PDBGeometry':
        from just_psf.parser.pdb import PDBParser
        return PDBParser(f).pdb()
