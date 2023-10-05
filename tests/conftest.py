import pathlib
import shutil
import pytest
import tempfile
import numpy

from just_psf.geometry import Geometry
from just_psf.topology import Topology

from tests import path_from_tests_files


@pytest.fixture(scope='module')
def tempdir():
    _tempdir = tempfile.mkdtemp()
    yield pathlib.Path(_tempdir)

    shutil.rmtree(_tempdir)


@pytest.fixture(scope='module')
def geometry_water():
    with path_from_tests_files(pathlib.Path('tests_files/H2O.xyz')).open() as f:
        geometry = Geometry.from_xyz(f)

    return geometry


@pytest.fixture(scope='module')
def geometry_fluoroethylene():
    with path_from_tests_files(pathlib.Path('tests_files/fluoroethylene.xyz')).open() as f:
        geometry = Geometry.from_xyz(f)

    return geometry


@pytest.fixture(scope='module')
def topology_water():
    """By hand"""

    return Topology(
        ['O', 'H', 'H'],
        ['O_3', 'H__HB', 'H__HB'],
        charges=[-.4, .2, .2],
        masses=[16.0, 1.01, 1.01],
        bonds=numpy.array([[1, 2], [1, 3]], dtype=int) - 1,
        angles=numpy.array([[2, 1, 3]], dtype=int) - 1,
        donors=numpy.array([[1, 2], [1, 3]], dtype=int) - 1,
        acceptors=numpy.array([[1, 0]], dtype=int) - 1
    )


@pytest.fixture(scope='module')
def topology_fluoroethylene():
    """By hand"""

    return Topology(
        ['F', 'C', 'C', 'H', 'H', 'H'],
        ['F_', 'C_2', 'C_2', 'H_', 'H_', 'H_'],
        masses=[19.0, 12.01, 12.01, 1.01, 1.01, 1.01],
        bonds=numpy.array([[1, 2], [2, 3], [2, 4], [3, 5], [3, 6]], dtype=int) - 1,
        angles=numpy.array([[1, 2, 3], [1, 2, 4], [2, 3, 5], [2, 3, 6], [4, 2, 3], [5, 3, 6]], dtype=int) - 1,
        dihedrals=numpy.array([[1, 2, 3, 5], [1, 2, 3, 6], [4, 2, 3, 5], [4, 2, 3, 6]], dtype=int) - 1,
        impropers=numpy.array([[2, 1, 4, 3], [3, 2, 5, 6]], dtype=int) - 1
    )


@pytest.fixture(scope='module')
def topology_fluoroethylene_psf():

    with path_from_tests_files(pathlib.Path('tests_files/fluoroethylene.psf')).open() as f:
        return Topology.from_psf(f)
