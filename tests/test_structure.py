import numpy
from io import StringIO

from just_psf.structure import Structure


def test_structure_water_ok(structure_water):
    assert len(structure_water) == 3


def test_structure_fluoroethylene_ok(structure_fluoroethylene):
    assert len(structure_fluoroethylene) == 6


def assert_structure_equals(topo1: Structure, topo2: Structure, restricted: bool = False):
    assert len(topo1) == len(topo2)

    assert topo1.atom_names == topo2.atom_names
    assert topo1.atom_types == topo2.atom_types

    if not restricted:
        assert topo1.seg_names == topo2.seg_names
        assert topo1.resi_ids == topo2.resi_ids
        assert topo1.resi_names == topo2.resi_names
        assert topo1.charges == topo2.charges
        assert topo1.masses == topo2.masses

    assert numpy.allclose(topo1.bonds, topo2.bonds)
    assert numpy.allclose(topo1.angles, topo2.angles)
    assert numpy.allclose(topo1.dihedrals, topo2.dihedrals)
    assert numpy.allclose(topo1.impropers, topo2.impropers)


def test_structure_read_psf_ok(structure_fluoroethylene_psf, structure_fluoroethylene):
    assert_structure_equals(structure_fluoroethylene, structure_fluoroethylene_psf, restricted=True)


def test_structure_write_psf_ok(structure_fluoroethylene_psf):
    f = StringIO()
    structure_fluoroethylene_psf.to_psf(f)

    f.seek(0)
    new_structure = Structure.from_psf(f)

    assert_structure_equals(structure_fluoroethylene_psf, new_structure)


def test_structure_write_psf_ext_ok(structure_fluoroethylene_psf):
    f = StringIO()
    structure_fluoroethylene_psf.to_psf(f, ['EXT'])

    f.seek(0)
    new_structure = Structure.from_psf(f)

    assert_structure_equals(structure_fluoroethylene_psf, new_structure)


def test_structure_write_psf_cheq_ok(structure_fluoroethylene_psf):
    f = StringIO()
    structure_fluoroethylene_psf.to_psf(f, ['CHEQ'])

    f.seek(0)
    new_structure = Structure.from_psf(f)

    assert_structure_equals(structure_fluoroethylene_psf, new_structure)


def test_structure_water_write_psf_ok(structure_water):
    assert len(structure_water) == 3

    f = StringIO()
    # ext xplor format required, because atom types are longer than 4 chars!
    structure_water.to_psf(f, ['EXT', 'XPLOR'])
    f.seek(0)

    new_structure = Structure.from_psf(f)

    assert structure_water.charges == new_structure.charges

    # check donors and acceptors:
    assert numpy.allclose(structure_water.bonds, new_structure.bonds)
    assert numpy.allclose(structure_water.angles, new_structure.angles)
    assert numpy.allclose(structure_water.donors, new_structure.donors)
    assert numpy.allclose(structure_water.acceptors, new_structure.acceptors)
