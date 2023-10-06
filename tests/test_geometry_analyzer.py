import numpy

from just_psf.geometry_analyzer import GeometryAnalyzer


def test_guess_bonds(geometry_fluoroethylene, structure_fluoroethylene):
    maker = GeometryAnalyzer(geometry_fluoroethylene)

    # check every bond was found
    for bond in structure_fluoroethylene.bonds:
        assert bond in maker.g.edges


def test_topology_make_water(geometry_water, structure_water):
    maker = GeometryAnalyzer(geometry_water)
    auto_topology = maker.structure()

    for bond in structure_water.bonds:
        assert bond in auto_topology.bonds or reversed(bond) in auto_topology.bonds

    for angle in structure_water.angles:
        assert angle in auto_topology.angles or reversed(angle) in auto_topology.angles


def test_topology_make_fluoroethylene(geometry_fluoroethylene, structure_fluoroethylene_psf):
    maker = GeometryAnalyzer(geometry_fluoroethylene)
    auto_topology = maker.structure()

    for bond in structure_fluoroethylene_psf.bonds:
        assert bond in auto_topology.bonds or reversed(bond) in auto_topology.bonds

    for angle in structure_fluoroethylene_psf.angles:
        assert angle in auto_topology.angles or reversed(angle) in auto_topology.angles

    for dihedral in structure_fluoroethylene_psf.dihedrals:
        assert dihedral in auto_topology.dihedrals or reversed(dihedral) in auto_topology.dihedrals


def test_topology_make_7waters(geometry_7waters):
    maker = GeometryAnalyzer(geometry_7waters)
    auto_topology = maker.structure()

    assert auto_topology.resi_ids == [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7]
    assert numpy.all(auto_topology.angles == [
        (1, 0, 2),
        (4, 3, 5),
        (7, 6, 8),
        (10, 9, 11),
        (13, 12, 14),
        (16, 15, 17),
        (19, 18, 20)
    ])
