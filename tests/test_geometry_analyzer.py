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
    auto_topology = maker.topology()

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
    auto_topology = maker.topology()

    assert auto_topology.masses.keys() == {'C', 'H', 'F'}
    assert auto_topology.autogenerate == {('ANGLE', 'DIHE')}

    assert len(auto_topology.residues) == 1
    assert auto_topology.residues[0].resi_name == 'MOL1'
    assert auto_topology.residues[0].atom_types == ['F', 'C', 'C', 'H', 'H', 'H']
    assert auto_topology.residues[0].atom_names == ['F1', 'C2', 'C3', 'H4', 'H5', 'H6']

    for bond in structure_fluoroethylene_psf.bonds:
        assert bond in auto_topology.residues[0].bonds or reversed(bond) in auto_topology.residues[0].bonds


def test_structure_make_7waters_ok(geometry_7waters):
    maker = GeometryAnalyzer(geometry_7waters)
    auto_structure = maker.structure()

    assert auto_structure.resi_ids == [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7]
    assert numpy.all(auto_structure.angles == [
        (1, 0, 2),
        (4, 3, 5),
        (7, 6, 8),
        (10, 9, 11),
        (13, 12, 14),
        (16, 15, 17),
        (19, 18, 20)
    ])
