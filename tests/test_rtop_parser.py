import io
import pytest


from just_psf.rtop_parser import RTopParser, RTopParseError


def test_parse_title_ok():
    title_lines = ['A title', 'on two lines']

    f = io.StringIO()
    f.write('*{}\n*\n'.format('\n*'.join(title_lines)))
    f.seek(0)

    parser = RTopParser(f)
    assert parser.title() == title_lines


def test_parse_wrong_title_ko():
    with pytest.raises(RTopParseError, match='expected `*`'):
        f = io.StringIO()
        f.write('*title\n')
        f.seek(0)

        RTopParser(f).title()


def test_parse_headers_ok():
    f = io.StringIO()

    # no title
    f.write('  19 1  \n')
    f.seek(0)

    parser = RTopParser(f)
    assert parser.title() == []
    assert parser.version() == (19, 1)

    # with title
    f.seek(0)
    f.write('*title\n*\n  36 1  \n')
    f.seek(0)

    parser = RTopParser(f)
    assert parser.title() == ['title']
    assert parser.version() == (36, 1)


def test_parse_decls():

    f = io.StringIO()

    masses = {'CH3': 12.01, 'HC': 1.01}
    f.write('\n'.join('MASS -1 {} {:.3f}'.format(k, v) for k, v in masses.items()))
    f.write('\n')

    decls = {'-O', '+N'}
    f.write('\n'.join('DECL {}'.format(d) for d in decls))
    f.write('\n')

    autos = {'ANGL', 'DIHE'}
    f.write('AUTO {}\n'.format(' '.join(autos)))

    defas = [('FIRS', 'NTER'), ('LAST', 'CTER')]
    f.write('\n'.join('DEFA {} {}'.format(a, b) for a, b in defas))
    f.write('\n')

    f.seek(0)

    parser = RTopParser(f)
    pmasses, pdecls, pdefas, pautos = parser.decls()

    assert pmasses == masses
    assert set(pdecls) == decls
    assert pdefas == defas
    assert set(pautos) == autos
