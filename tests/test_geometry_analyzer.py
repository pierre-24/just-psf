import numpy

from just_psf.geometry_analyzer import GeometryAnalyzer


def test_guess_bonds_ok(geometry_fluoroethylene, structure_fluoroethylene):
    maker = GeometryAnalyzer(geometry_fluoroethylene)

    # check every bond was found
    for bond in structure_fluoroethylene.bonds:
        assert bond in maker.g.edges


def test_structure_make_water_ok(geometry_water, structure_water):
    maker = GeometryAnalyzer(geometry_water)
    auto_structure = maker.structure()

    for bond in structure_water.bonds:
        assert bond in auto_structure.bonds or reversed(bond) in auto_structure.bonds

    for angle in structure_water.angles:
        assert angle in auto_structure.angles or reversed(angle) in auto_structure.angles


def test_topology_make_water(geometry_water, structure_water):
    maker = GeometryAnalyzer(geometry_water)
    auto_topology = maker.topologies()

    assert auto_topology.masses.keys() == {'H', 'O'}
    assert auto_topology.autogenerate == {('ANGLE', 'DIHE')}

    assert len(auto_topology.residues) == 1
    assert auto_topology.residues[0].resi_name == 'MOL1'
    assert auto_topology.residues[0].atom_types == ['O', 'H', 'H']
    assert auto_topology.residues[0].atom_names == ['O1', 'H2', 'H3']

    for bond in structure_water.bonds:
        assert bond in auto_topology.residues[0].bonds or reversed(bond) in auto_topology.residues[0].bonds


def test_structure_make_fluoroethylene_ok(geometry_fluoroethylene, structure_fluoroethylene_psf):
    maker = GeometryAnalyzer(geometry_fluoroethylene)
    auto_structure = maker.structure()

    for bond in structure_fluoroethylene_psf.bonds:
        assert bond in auto_structure.bonds or reversed(bond) in auto_structure.bonds

    for angle in structure_fluoroethylene_psf.angles:
        assert angle in auto_structure.angles or reversed(angle) in auto_structure.angles

    for dihedral in structure_fluoroethylene_psf.dihedrals:
        assert dihedral in auto_structure.dihedrals or reversed(dihedral) in auto_structure.dihedrals


def test_topology_make_fluoroethylene_ok(geometry_fluoroethylene, structure_fluoroethylene_psf):
    maker = GeometryAnalyzer(geometry_fluoroethylene)
    auto_topology = maker.topologies()

    assert auto_topology.masses.keys() == {'C', 'H', 'F'}
    assert auto_topology.autogenerate == {('ANGLE', 'DIHE')}

    assert len(auto_topology.residues) == 1
    assert auto_topology.residues[0].resi_name == 'MOL1'
    assert auto_topology.residues[0].atom_types == ['F', 'C', 'C', 'H', 'H', 'H']
    assert auto_topology.residues[0].atom_names == ['F1', 'C2', 'C3', 'H4', 'H5', 'H6']

    for bond in structure_fluoroethylene_psf.bonds:
        assert bond in auto_topology.residues[0].bonds or reversed(bond) in auto_topology.residues[0].bonds


def test_structure_make_7waters_ok(geometry_7waters, structure_7water_psf):
    maker = GeometryAnalyzer(geometry_7waters)
    auto_structure = maker.structure()

    assert auto_structure.resi_ids == structure_7water_psf.resi_ids
    assert numpy.allclose(auto_structure.bonds, structure_7water_psf.bonds)
    assert numpy.allclose(auto_structure.angles, structure_7water_psf.angles)
