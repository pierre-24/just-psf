from typing import List, Dict, Set, Tuple, Optional
from numpy.typing import NDArray


class ResidueTopology:
    """
    Single residue topology.
    Note that `bonds` contains a tuple of two indices, which refers to the corresponding atom if >=0, but refer to
    the `declarations` if <0.
    """
    def __init__(
        self,
        resi_name: str,
        atom_names: List[str],
        atom_types: List[str],
        atom_charges: List[float],
        bonds: NDArray[int],
        resi_charge: float = .0,
    ):
        assert len(atom_types) == len(atom_names)
        assert len(atom_charges) == len(atom_charges)
        assert bonds.shape[1] == 2

        self.resi_name = resi_name
        self.resi_charge = resi_charge
        self.atom_names = atom_names
        self.atom_types = atom_types
        self.atom_charges = atom_charges
        self.bonds = bonds

    def __len__(self):
        return len(self.atom_names)

    def as_rtop(self, declarations: List[str]) -> str:
        # resi
        r = 'RESI {:4} {: .2f}\nGROUP\n'.format(self.resi_name, self.resi_charge)

        # atoms
        for i in range(len(self)):
            r += 'ATOM {:4} {:4} {: .2f}\n'.format(self.atom_names[i], self.atom_types[i], self.atom_charges[i])

        # bonds
        r += '!'
        for i in range(self.bonds.shape[0]):
            if i % 4 == 0:
                r += '\nBOND'

            a1, a2 = self.bonds[i]
            na1 = self.atom_names[a1] if a1 >= 0 else declarations[-a1 - 1]
            na2 = self.atom_names[a2] if a2 >= 0 else declarations[-a2 - 1]
            r += ' {:4} {:4}'.format(na1, na2)

        r += '\n'

        return r


class ResidueTopologies:
    def __init__(
        self,
        masses: Dict[str, float],
        defaults: Optional[Set[Tuple[str, str]]] = None,
        autogenerate: Optional[Set[tuple]] = None,
        declarations: Optional[List[str]] = None,
        residues: Optional[List[ResidueTopology]] = None,
    ):
        self.residues = residues if residues is not None else []
        self.masses = masses
        self.defaults = defaults if defaults is not None else set()
        self.autogenerate = autogenerate if autogenerate is not None else set()
        self.declarations = declarations if declarations is not None else []

    def as_rtop(self, version: int = 19) -> str:
        r = '* Topology generated by `just_psf.residue_topology.ResidueTopologies.as_rtop()`\n*\n19 1\n\n'

        # masses
        for name, mass in self.masses.items():
            r += 'MASS -1 {:4} {:7.3f}\n'.format(name, mass)

        r += '\n'

        # defaults
        if len(self.defaults) > 0:
            r += 'DEFA'
            for default in self.defaults:
                r += ' {} {}'.format(*default)
            r += '\n'

        # autogenerate
        for autogen in self.autogenerate:
            r += 'AUTO {}\n'.format(' '.join(autogen))

        # decls
        for name in self.declarations:
            r += 'DECL {}\n'.format(name)

        r += '\n'

        for residue in self.residues:
            r += residue.as_rtop(self.declarations)
            r += '\n'

        # the end:
        r += 'END\n'

        return r
