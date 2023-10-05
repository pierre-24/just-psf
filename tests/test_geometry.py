import pathlib

import numpy

from tests import path_from_tests_files


def test_water_ok(geometry_water):
    with path_from_tests_files(pathlib.Path('tests_files/H2O.xyz')).open() as f:
        content = f.read()

        lines = content.splitlines()

    assert len(geometry_water) == int(content[0])
    assert geometry_water.symbols == ['O', 'H', 'H']
    assert numpy.allclose(geometry_water.positions[0], [float(x) for x in lines[2].split()[1:]])
    assert numpy.allclose(geometry_water.positions[1], [float(x) for x in lines[3].split()[1:]])
    assert numpy.allclose(geometry_water.positions[2], [float(x) for x in lines[4].split()[1:]])
