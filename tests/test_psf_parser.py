from io import StringIO
import random
from typing import TextIO

import numpy
import pytest

from just_psf.parser.psf import PSFParser, PSFParseError


def random_indices(f: TextIO, intsize: int, n: int, indices_per: int, elements_per: int, excess: int = 0):
    seq = []
    intfmt = '{{:{}}}'.format(intsize)
    for i in range((n + excess) * indices_per):
        if i != 0 and i % (elements_per * indices_per) == 0:
            f.write('\n')

        num = random.randint(0, 255)
        f.write(intfmt.format(num + 1))  # indices are 1-based in PSF
        seq.append(num)

    f.write('\n')

    return numpy.int64(seq).reshape(n + excess, indices_per)


def test_read_indices_ok():

    def assert_seq_eq(intsize: int, n: int, indices_per: int, elements_per: int):
        f = StringIO()
        seq = random_indices(f, intsize, n, indices_per, elements_per)
        f.seek(0)

        assert numpy.allclose(PSFParser(f).parse_indices(intsize, n, indices_per, elements_per), seq)

    # "bond" section (elements of 2 indices each, 4 elements per line) of PSF:
    assert_seq_eq(8, 6, 2, 4)

    # "angle" section (elements of 3 indices each, 3 elements per line) of EXT PSF:
    assert_seq_eq(10, 8, 3, 3)


def test_read_bogus_indices_ko():

    def assert_raises(intsize: int, n: int, indices_per: int, elements_per: int, excess: int, match: str):
        f = StringIO()
        random_indices(f, intsize, n, indices_per, elements_per, excess)
        f.seek(0)

        with pytest.raises(PSFParseError, match=match):
            PSFParser(f).parse_indices(intsize, n, indices_per, elements_per)

    # too short
    assert_raises(8, 6, 2, 4, -1, 'incorrect number of indices')
    assert_raises(8, 10, 2, 4, -6, 'not enough')  # exactly one line

    # too long
    assert_raises(8, 6, 2, 4, 1, 'incorrect number of indices')
    assert_raises(8, 6, 2, 4, 2, 'incorrect number of indices')
    assert_raises(8, 8, 2, 4, 2, 'too much data')  # exactly one line


def random_atoms(f: TextIO, n: int, fmt: str = '{:>8d} {:4} {:4d} {:4} {:4} {:4} {:>14.6f}{:>14.6f}{:8d}'):
    atoms = []
    for i in range(n):
        atoms.append((
            i + 1,
            'SYS',
            1,
            'X',
            random.choice(['H', 'C', 'N', 'O', 'F']),
            'X',
            random.randrange(-1, 1),
            random.randrange(0, 20),
            0,
        ))
        f.write(fmt.format(*atoms[-1]))
        f.write('\n')

    return atoms


def test_read_atoms_ok():
    N = 5

    f = StringIO()
    atoms = random_atoms(f, N)
    f.seek(0)

    fid, segn, resi, resn, anam, atyp, chrg, mass, fixd = PSFParser(f).parse_atoms(N, [])

    assert fid == 1
    assert segn == list(a[1] for a in atoms)
    assert resi == list(a[2] for a in atoms)
    assert resn == list(a[3] for a in atoms)
    assert anam == list(a[4] for a in atoms)
    assert atyp == list(a[5] for a in atoms)
    assert chrg == list(a[6] for a in atoms)
    assert mass == list(a[7] for a in atoms)
    assert fixd == list(a[8] for a in atoms)


def test_read_too_much_atoms_ko():
    N = 5

    f = StringIO()
    random_atoms(f, N)
    f.seek(0)

    with pytest.raises(PSFParseError, match='too much'):
        PSFParser(f).parse_atoms(N - 1, [])


def test_read_not_enough_atoms_ko():
    N = 5

    f = StringIO()
    random_atoms(f, N)
    f.seek(0)

    with pytest.raises(PSFParseError, match='not enough'):
        PSFParser(f).parse_atoms(N + 1, [])


def test_parse_incorrect_psf_ko():
    def parse(inp: str):
        f = StringIO()
        f.write(inp)
        f.seek(0)
        PSFParser(f).structure()

    # missing flags line
    with pytest.raises(PSFParseError, match='not a PSF'):
        parse('       0 !NBOND\n')

    # missing PSF flag
    with pytest.raises(PSFParseError, match='not a PSF'):
        parse('EXT\n       0 !NBOND\n')

    # expected empty
    with pytest.raises(PSFParseError, match='expected empty'):
        parse('PSF\n       0 !NBOND\n')

    # too much empty
    with pytest.raises(PSFParseError, match='expected section'):
        parse('PSF\n\n\n       0 !NBOND\n')

    # not the correct integer format
    with pytest.raises(PSFParseError, match='incorrectly formatted'):
        parse('PSF\n\n           3 !NBOND: bonds')

    # missing NATOM
    with pytest.raises(PSFParseError, match='missing NATOM'):
        parse('PSF\n\n       0 !NBOND\n')
