from typing import List, Dict, Set, Tuple
from numpy.typing import NDArray


class ResidueTopology:
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


class ResidueTopologies:
    def __init__(
        self,
        masses: Dict[str, float],
        defaults: Set[Tuple[str, str]],
        autogenerate: Set[str],
        declarations: Set[str],
        residues: List[ResidueTopology]
    ):
        self.residues = residues
        self.masses = masses
        self.defaults = defaults
        self.autogenerate = autogenerate
        self.declarations = declarations
