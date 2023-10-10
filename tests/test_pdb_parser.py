import io
import numpy

from just_psf.geometry import PDBGeometry


def test_parse_pdb_ok(geometry_7waters_pdb, structure_7water_psf, geometry_7waters):
    assert geometry_7waters_pdb.symbols == geometry_7waters.symbols
    assert numpy.allclose(geometry_7waters_pdb.positions, geometry_7waters.positions, atol=1e-3)

    assert geometry_7waters_pdb.resi_ids == structure_7water_psf.resi_ids
    assert geometry_7waters_pdb.atom_names == structure_7water_psf.atom_names
    assert geometry_7waters_pdb.resi_names == ['HOH'] * 21
    assert geometry_7waters_pdb.seg_names == [''] * 21


def test_write_pdb_ok(geometry_7waters_pdb):
    f = io.StringIO()

    geometry_7waters_pdb.to_pdb(f)
    f.seek(0)

    geom = PDBGeometry.from_pdb(f)

    assert geometry_7waters_pdb.resi_ids == geom.resi_ids
    assert geometry_7waters_pdb.resi_names == geom.resi_names
    assert geometry_7waters_pdb.atom_names == geom.atom_names
    assert geometry_7waters_pdb.seg_names == geom.seg_names
