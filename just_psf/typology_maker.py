from typing import Union

import networkx
import numpy
from scipy.spatial import distance_matrix
from typing import Iterable
import queue

from just_psf import logger
from just_psf.geometry import Geometry
from just_psf.topology import Topology


l_logger = logger.getChild(__name__)


COVALENT_RADII = {
    # from 10.1039/B801115J
    'H': 0.31,
    'He': 0.28,
    'Li': 1.28,
    'Be': 0.96,
    'B': 0.84,
    'C': 0.76,
    'N': 0.71,
    'O': 0.66,
    'F': 0.57,
    'Ne': 0.58,
    'Na': 1.66,
    'Mg': 1.41,
    'Al': 1.21,
    'Si': 1.11,
    'P': 1.07,
    'S': 1.05,
    'Cl': 1.02,
    'Ar': 1.06,
    'K': 2.03,
    'Ca': 1.76,
    'Sc': 1.70,
    'Ti': 1.60,
    'V': 1.53,
    'Cr': 1.39,
    'Mn': 1.39,
    'Fe': 1.32,
    'Co': 1.26,
    'Ni': 1.24,
    'Cu': 1.32,
    'Zn': 1.22,
    'Ga': 1.22,
    'Ge': 1.20,
    'As': 1.19,
    'Se': 1.20,
    'Br': 1.20,
    'Kr': 1.16,
    'Rb': 2.20,
    'Sr': 1.95,
    'Y': 1.90,
    'Zr': 1.75,
    'Nb': 1.64,
    'Mo': 1.54,
    'Tc': 1.47,
    'Ru': 1.46,
    'Rh': 1.42,
    'Pd': 1.39,
    'Ag': 1.45,
    'Cd': 1.44,
    'In': 1.42,
    'Sn': 1.39,
    'Sb': 1.39,
    'Te': 1.38,
    'I': 1.39,
    'Xe': 1.40,
    'Cs': 2.44,
    'Ba': 2.15,
    'La': 2.07,
    'Ce': 2.04,
    'Pr': 2.03,
    'Nd': 2.01,
    'Pm': 1.99,
    'Sm': 1.98,
    'Eu': 1.98,
    'Gd': 1.96,
    'Tb': 1.94,
    'Dy': 1.92,
    'Ho': 1.92,
    'Er': 1.89,
    'Tm': 1.90,
    'Yb': 1.87,
    'Lu': 1.87,
    'Hf': 1.75,
    'Ta': 1.70,
    'W': 1.62,
    'Re': 1.51,
    'Os': 1.44,
    'Ir': 1.41,
    'Pt': 1.36,
    'Au': 1.36,
    'Hg': 1.32,
    'Tl': 1.45,
    'Pb': 1.46,
    'Bi': 1.48,
    'Po': 1.40,
    'At': 1.50,
    'Rn': 1.50,
    'Fr': 2.60,
    'Ra': 2.21,
    'Ac': 2.15,
    'Th': 2.06,
    'Pa': 2.00,
    'U': 1.96,
    'Np': 1.90,
    'Pu': 1.87,
    'Am': 1.80,
    'Cm': 1.69,
}


def find_subgraphs(g: networkx.Graph, k: int) -> Iterable[tuple]:
    """Find all subgraphs of `g` of size `k` that does not form nor contain a loop.
    Also report only one of the two symmetric paths (i.e., `i->j->k` and not `k->j->i`).
    """

    def paths(root: int) -> Iterable[tuple]:
        """Yield all path of size `k` from `root` in `g` that does not form nor contain a loop.
        """

        # breadth-first search:
        q = queue.SimpleQueue()
        q.put((root,))

        while not q.empty():
            p = q.get()
            if len(p) < k:
                for j in g.adj[p[-1]]:
                    if j in p:  # avoid loops
                        continue

                    q.put(p + (j,))
            else:
                yield p

    for i in g.nodes:
        for subgraph in paths(i):
            if subgraph[-1] < i:  # avoid reversed
                continue

            yield subgraph


class TopologyMaker:
    def __init__(self, geometry: Union[str, Geometry], threshold: float = 1.1):
        if type(geometry) is str:
            with open(geometry) as f:
                self.geometry = Geometry.from_xyz(f)
        elif type(geometry) is Geometry:
            self.geometry = geometry
        else:
            raise TypeError('geometry')

        self.g = networkx.Graph()
        self.g.add_nodes_from(range(len(geometry)))

        # get bonds
        self._guess_bonds(threshold)

    def _guess_bonds(self, threshold: float = 1.1):
        """Guess which atom are linked to which using a distance matrix.
        May lead to incorrect results for strange bonds (e.g., metalic)
        """
        l_logger.info('compute distances')
        distances = distance_matrix(self.geometry.positions, self.geometry.positions)
        l_logger.info('assign bonds')
        for i in range(len(self.geometry)):
            cri = COVALENT_RADII[self.geometry.symbols[i]]
            for j in range(i + 1, len(self.geometry)):
                crj = COVALENT_RADII[self.geometry.symbols[j]]
                if distances[i, j] < threshold * (cri + crj):
                    self.g.add_edge(i, j)

    def topology(self) -> Topology:
        """
        Get the corresponding topology
        """

        angles = []
        dihedrals = []

        mol_ids = [0] * len(self.geometry)

        i = 0
        for indices in networkx.connected_components(self.g):
            i += 1
            l_logger.info('* check molecule {} (size = {} atoms)'.format(i, len(indices)))
            for index in indices:
                mol_ids[index] = i

            subgraph = self.g.subgraph(indices)
            l_logger.info('get angles')
            if len(subgraph.nodes) > 2:
                angles.extend(find_subgraphs(subgraph, 3))
            l_logger.info('get dihedrals')
            if len(subgraph.nodes) > 4:
                dihedrals.extend(find_subgraphs(subgraph, 4))

        return Topology(
            atom_names=self.geometry.symbols,
            atom_types=self.geometry.symbols,
            atom_ids=list(range(len(self.geometry))),
            resi_ids=mol_ids,
            bonds=numpy.array(self.g.edges),
            angles=numpy.array(angles) if len(angles) > 0 else None,
            dihedrals=numpy.array(dihedrals) if len(dihedrals) > 0 else None
        )
