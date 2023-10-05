import numpy
from io import StringIO

from just_psf.topology import Topology


def test_topology_water_ok(topology_water):
    assert len(topology_water) == 3


def test_topology_fluoroethylene_ok(topology_fluoroethylene):
    assert len(topology_fluoroethylene) == 6


def assert_fluoroethylene_topologies_equals(topo1: Topology, topo2: Topology, restricted: bool = False):
    assert len(topo1) == len(topo2)

    assert topo1.atom_names == topo2.atom_names
    assert topo1.atom_types == topo2.atom_types

    if not restricted:
        assert topo1.atom_ids == topo2.atom_ids
        assert topo1.seg_names == topo2.seg_names
        assert topo1.resi_ids == topo2.resi_ids
        assert topo1.resi_names == topo2.resi_names
        assert topo1.charges == topo2.charges
        assert topo1.masses == topo2.masses

    assert numpy.allclose(topo1.bonds, topo2.bonds)
    assert numpy.allclose(topo1.angles, topo2.angles)
    assert numpy.allclose(topo1.dihedrals, topo2.dihedrals)
    assert numpy.allclose(topo1.impropers, topo2.impropers)


def test_topology_read_psf_ok(topology_fluoroethylene_psf, topology_fluoroethylene):
    assert_fluoroethylene_topologies_equals(topology_fluoroethylene, topology_fluoroethylene_psf, restricted=True)


def test_topology_write_psf_ok(topology_fluoroethylene_psf):
    f = StringIO()
    topology_fluoroethylene_psf.to_psf(f)

    f.seek(0)
    new_topology = Topology.from_psf(f)

    assert_fluoroethylene_topologies_equals(topology_fluoroethylene_psf, new_topology)


def test_topology_write_psf_ext_ok(topology_fluoroethylene_psf):
    f = StringIO()
    topology_fluoroethylene_psf.to_psf(f, ['EXT'])

    f.seek(0)
    new_topology = Topology.from_psf(f)

    assert_fluoroethylene_topologies_equals(topology_fluoroethylene_psf, new_topology)


def test_topology_write_psf_cheq_ok(topology_fluoroethylene_psf):
    f = StringIO()
    topology_fluoroethylene_psf.to_psf(f, ['CHEQ'])

    f.seek(0)
    new_topology = Topology.from_psf(f)

    assert_fluoroethylene_topologies_equals(topology_fluoroethylene_psf, new_topology)


def test_topology_water_write_psf_ok(topology_water):
    assert len(topology_water) == 3

    f = StringIO()
    topology_water.to_psf(f, ['EXT', 'XPLOR'])  # ext xplor format required, because atom types are longer than 4 chars!
    f.seek(0)
    new_topology = Topology.from_psf(f)

    assert topology_water.charges == new_topology.charges

    # check donors and acceptors:
    assert numpy.allclose(topology_water.bonds, new_topology.bonds)
    assert numpy.allclose(topology_water.angles, new_topology.angles)
    assert numpy.allclose(topology_water.donors, new_topology.donors)
    assert numpy.allclose(topology_water.acceptors, new_topology.acceptors)
