from typing import Union, List, Tuple

import networkx
import numpy
from scipy.spatial import distance_matrix
from typing import Iterable
import queue

from just_psf import logger
from just_psf.geometry import Geometry
from just_psf.residue_topology import Topologies, ResidueTopology
from just_psf.structure import Structure


l_logger = logger.getChild(__name__)


COVALENT_RADII = {  # from 10.1039/B801115J
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

ATOMIC_WEIGHTS = {  # from https://iupac.qmul.ac.uk/AtWt/
    'H': 1.008,
    'He': 4.003,
    'Li': 6.940,
    'Be': 9.012,
    'B': 10.810,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    'F': 18.998,
    'Ne': 20.180,
    'Na': 22.990,
    'Mg': 24.305,
    'Al': 26.982,
    'Si': 28.085,
    'P': 30.974,
    'S': 32.060,
    'Cl': 35.450,
    'Ar': 39.950,
    'K': 39.098,
    'Ca': 40.078,
    'Sc': 44.956,
    'Ti': 47.867,
    'V': 50.942,
    'Cr': 51.996,
    'Mn': 54.938,
    'Fe': 55.845,
    'Co': 58.933,
    'Ni': 58.693,
    'Cu': 63.546,
    'Zn': 65.380,
    'Ga': 69.723,
    'Ge': 72.630,
    'As': 74.922,
    'Se': 78.971,
    'Br': 79.904,
    'Kr': 83.798,
    'Rb': 85.468,
    'Sr': 87.620,
    'Y': 88.906,
    'Zr': 91.224,
    'Nb': 92.906,
    'Mo': 95.950,
    'Tc': 97.000,
    'Ru': 101.070,
    'Rh': 102.905,
    'Pd': 106.420,
    'Ag': 107.868,
    'Cd': 112.414,
    'In': 114.818,
    'Sn': 118.710,
    'Sb': 121.760,
    'Te': 127.600,
    'I': 126.904,
    'Xe': 131.293,
    'Cs': 132.905,
    'Ba': 137.327,
    'La': 138.905,
    'Ce': 140.116,
    'Pr': 140.908,
    'Nd': 144.242,
    'Pm': 145.000,
    'Sm': 150.360,
    'Eu': 151.964,
    'Gd': 157.250,
    'Tb': 158.925,
    'Dy': 162.500,
    'Ho': 164.930,
    'Er': 167.259,
    'Tm': 168.934,
    'Yb': 173.045,
    'Lu': 174.967,
    'Hf': 178.486,
    'Ta': 180.948,
    'W': 183.840,
    'Re': 186.207,
    'Os': 190.230,
    'Ir': 192.217,
    'Pt': 195.084,
    'Au': 196.967,
    'Hg': 200.592,
    'Tl': 204.380,
    'Pb': 207.200,
    'Bi': 208.980,
    'Po': 209.000,
    'At': 210.000,
    'Rn': 222.000,
    'Fr': 223.000,
    'Ra': 226.000,
    'Ac': 227.000,
    'Th': 232.038,
    'Pa': 231.036,
    'U': 238.029,
    'Np': 237.000,
    'Pu': 244.000,
    'Am': 243.000,
    'Cm': 247.000,
    'Bk': 247.000,
    'Cf': 251.000,
    'Es': 252.000,
    'Fm': 257.000,
    'Md': 258.000,
    'No': 259.000,
    'Lr': 262.000,
    'Rf': 267.000,
    'Db': 270.000,
    'Sg': 269.000,
    'Bh': 270.000,
    'Hs': 270.000,
    'Mt': 278.000,
    'Ds': 281.000,
    'Rg': 281.000,
    'Cn': 285.000,
    'Nh': 286.000,
    'Fl': 289.000,
    'Mc': 289.000,
    'Lv': 293.000,
    'Ts': 293.000,
    'Og': 294.000,
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


class MolecularSubgraph:
    """A subgraph which represent a "molecule", i.e., a connected component in said graph.
    """

    def __init__(self, subgraph: networkx.Graph):
        self.subgraph = subgraph

    def autogenerate_angles_dihedrals(self) -> Tuple[List[Tuple[int, int, int]], List[Tuple[int, int, int, int]]]:
        l_logger.debug('Generate angles and dihedrals')

        angles = []
        dihedrals = []

        if len(self.subgraph.nodes) > 2:
            angles = list(find_subgraphs(self.subgraph, 3))

        if len(self.subgraph.nodes) > 3:
            dihedrals = list(find_subgraphs(self.subgraph, 4))

        return angles, dihedrals

    def __len__(self):
        return len(self.subgraph.nodes)


class GeometryAnalyzer:
    def __init__(self, geometry: Union[str, Geometry], threshold: float = 1.1):
        if type(geometry) is str:
            with open(geometry) as f:
                self.geometry = Geometry.from_xyz(f)
        elif type(geometry) is Geometry:
            self.geometry = geometry
        else:
            raise TypeError('geometry')

        # create graph
        self.g = networkx.Graph()
        self.g.add_nodes_from((i, {'symbol': self.geometry.symbols[i]}) for i in range(len(geometry)))
        self._guess_bonds(threshold)

        # get and analyze connected components
        self.atom_names = [''] * len(self.geometry)
        self.resi_ids = [0] * len(self.geometry)

        self.resi_isomorphic_to = {}
        self.atom_isomorphic_to = [-1] * len(self.geometry)

        self.uniq_residues = []
        current_resi_id = -1
        for indices in networkx.connected_components(self.g):
            current_resi_id += 1
            subgraph = self.g.subgraph(indices)

            # tries to match the current residue to another
            uniq_resi_id = -1
            for i, uniq in enumerate(self.uniq_residues):
                gm = networkx.isomorphism.GraphMatcher(
                    uniq.subgraph,
                    subgraph,
                    node_match=networkx.isomorphism.categorical_node_match('symbol', 'X')
                )
                if gm.is_isomorphic():
                    uniq_resi_id = i
                    mapping = gm.mapping
                    break

            if uniq_resi_id < 0:
                uniq_resi_id = len(self.uniq_residues)
                mapping = dict((i, i) for i in indices)
                self.uniq_residues.append(MolecularSubgraph(subgraph))

            # store isomorphism
            if uniq_resi_id not in self.resi_isomorphic_to:
                self.resi_isomorphic_to[uniq_resi_id] = []

            self.resi_isomorphic_to[uniq_resi_id].append(mapping)

            # fill resi_id and atom_names
            for uresi_ai, resi_ai in mapping.items():
                self.atom_isomorphic_to[resi_ai] = uresi_ai
                self.resi_ids[resi_ai] = current_resi_id + 1
                self.atom_names[resi_ai] = '{}{}'.format(
                    self.geometry.symbols[uresi_ai],
                    uresi_ai + 1
                )

        l_logger.info('found {} connected component(s) and {} unique(s)'.format(
            current_resi_id, len(self.uniq_residues)))

    def _guess_bonds(self, threshold: float = 1.1):
        """
        Guess which atom are linked to which using a distance matrix.
        May lead to incorrect results for strange bonds (e.g., metalic)
        """
        l_logger.debug('compute distances')
        distances = distance_matrix(self.geometry.positions, self.geometry.positions)

        l_logger.debug('assign bonds')
        for i in range(len(self.geometry)):
            cri = COVALENT_RADII[self.geometry.symbols[i]]
            for j in range(i + 1, len(self.geometry)):
                crj = COVALENT_RADII[self.geometry.symbols[j]]
                if distances[i, j] < threshold * (cri + crj):
                    self.g.add_edge(i, j)

    def structure(self, seg_name: str = 'SYS') -> Structure:
        """
        Get the corresponding structure.
        Each connected component is considered as a residue.
        """

        angles = []
        dihedrals = []

        resi_names = ['X'] * len(self.geometry)

        for i, component in enumerate(self.uniq_residues):
            resi_angs, resi_dihe = component.autogenerate_angles_dihedrals()
            for mp in self.resi_isomorphic_to[i]:
                # assign resi name
                for ai in mp.values():
                    resi_names[ai] = 'RES{}'.format(i + 1)

                # add angles using mapping
                angles.extend((mp[i], mp[j], mp[k]) for i, j, k in resi_angs)
                dihedrals.extend((mp[i], mp[j], mp[k], mp[l]) for i, j, k, l in resi_dihe)

        return Structure(
            seg_names=[seg_name] * len(self.geometry),
            atom_types=self.geometry.symbols,
            atom_names=self.atom_names,
            resi_ids=self.resi_ids,
            resi_names=resi_names,
            masses=[ATOMIC_WEIGHTS[s] for s in self.geometry.symbols],
            bonds=numpy.array(self.g.edges),
            angles=numpy.array(angles) if len(angles) > 0 else None,
            dihedrals=numpy.array(dihedrals) if len(dihedrals) > 0 else None
        )

    def topologies(self) -> Topologies:
        """
        Get a set of topologies.
        Each independent component is converted in a residue.
        """

        uniq_elements = set(self.geometry.symbols)

        residues = []
        i = 0
        for residue in self.uniq_residues:
            i += 1
            residues.append(ResidueTopology(
                resi_name='RES{}'.format(i),
                resi_charge=.0,
                atom_types=[self.geometry.symbols[i] for i in residue.subgraph.nodes],
                atom_names=[self.atom_names[i] for i in residue.subgraph.nodes],
                atom_charges=[.0] * len(residue),
                bonds=numpy.array(residue.subgraph.edges)
            ))

        return Topologies(
            masses=dict((k, ATOMIC_WEIGHTS[k]) for k in uniq_elements),
            autogenerate={('ANGLE', 'DIHE')},
            defaults={('FIRST', 'NONE'), ('LAST', 'NONE')},
            residues=residues,
        )
