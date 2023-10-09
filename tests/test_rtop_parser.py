import io
import pathlib

import numpy
import pytest


from just_psf.parser.rtop import RTopParser, RTopParseError
from tests import path_from_tests_files


def test_parse_title_ok():
    title_lines = ['A title', 'on two lines']
    parser = RTopParser('*{}\n*\n'.format('\n*'.join(title_lines)))
    assert parser.title() == title_lines


def test_parse_wrong_title_ko():
    with pytest.raises(RTopParseError, match='expected `*`'):
        RTopParser('*title\n').title()


def test_parse_headers_ok():
    # no title
    parser = RTopParser('  19 1  \n')
    assert parser.title() == []
    assert parser.version() == (19, 1)

    # with title
    parser = RTopParser('*title\n*\n  36 1  \n')
    assert parser.title() == ['title']
    assert parser.version() == (36, 1)

    # with everything else
    f = io.StringIO()
    f.write('* Title\n*\n19 1\n')

    masses = {'CH3': 12.01, 'HC': 1.01}
    f.write('\n'.join('MASS -1 {} {:.3f}'.format(k, v) for k, v in masses.items()))
    f.write('\n')

    decls = ['-C', '+C']
    f.write('\n'.join('DECL {}'.format(d) for d in decls))
    f.write('\n')

    autos = {('ANGL', 'DIHE'), ('ANGLE', 'DIHE', 'PATCH')}
    f.write('\n'.join('AUTO ' + ' '.join(a for a in x) for x in autos))
    f.write('\n')

    defas = {('FIRS', 'NTER'), ('LAST', 'CTER')}
    f.write('\n'.join('DEFA {} {}'.format(a, b) for a, b in defas))
    f.write('\n')

    f.write('END')

    f.seek(0)

    parser = RTopParser(f)
    topology = parser.topologies()

    assert topology.declarations == decls
    assert topology.defaults == defas
    assert topology.masses == masses
    assert topology.autogenerate == autos


def assert_residue_equals(residue, residue2):
    assert residue.resi_name == residue2.resi_name
    assert residue.resi_charge == residue2.resi_charge
    assert residue.atom_names == residue2.atom_names
    assert residue.atom_types == residue2.atom_types
    assert residue.atom_charges == residue2.atom_charges
    assert numpy.allclose(residue.bonds, residue2.bonds)


def test_parse_residue_ok():
    allowed_types = {'CH3', 'HC'}
    declarations = ['-C', '+C']
    parser = RTopParser('RESI TEST 0.0\nATOM C CH3 0.0\nATOM H1 HC 0.0\nATOM H2 HC 0.0\nBOND C H1 C H2 -C C C +C\nEND')
    residue = parser.residue(allowed_types, declarations)

    assert residue.resi_name == 'TEST'
    assert residue.resi_charge == .0
    assert residue.atom_names == ['C', 'H1', 'H2']
    assert residue.atom_types == ['CH3', 'HC', 'HC']
    assert residue.atom_charges == [.0, .0, .0]
    assert numpy.allclose(residue.bonds, [(0, 1), (0, 2), (-1, 0), (0, -2)])

    assert parser.current_token.value == 'END'

    rtop = residue.as_rtop(declarations)
    parser = RTopParser(rtop)

    residue2 = parser.residue(allowed_types, declarations)
    assert_residue_equals(residue, residue2)


def test_parse_topology_ok():
    with path_from_tests_files(pathlib.Path('tests_files/topology.tpr')).open() as f:
        parser = RTopParser(f)
        topology = parser.topologies()

    assert [r.resi_name for r in topology.residues] == ['ALA', 'ARG']

    parser = RTopParser(topology.as_rtop())
    topology2 = parser.topologies()

    assert topology.masses.keys() == topology2.masses.keys()
    assert topology.defaults == topology2.defaults
    assert topology.declarations == topology2.declarations
    assert topology.autogenerate == topology2.autogenerate

    assert_residue_equals(topology.residues[0], topology2.residues[0])
    assert_residue_equals(topology.residues[1], topology2.residues[1])
