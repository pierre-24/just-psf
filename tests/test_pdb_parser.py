import numpy


def test_parse_pdb(geometry_7waters_pdb, structure_7water_psf, geometry_7waters):
    assert geometry_7waters_pdb.symbols == geometry_7waters.symbols
    assert numpy.allclose(geometry_7waters_pdb.positions, geometry_7waters.positions, atol=1e-3)

    assert geometry_7waters_pdb.resi_ids == structure_7water_psf.resi_ids
    assert geometry_7waters_pdb.atom_types == structure_7water_psf.atom_types
    assert geometry_7waters_pdb.resi_names == ['HOH'] * 21
    assert geometry_7waters_pdb.seg_names == [''] * 21
