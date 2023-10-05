from just_psf.typology_maker import TopologyMaker


def test_guess_bonds(geometry_fluoroethylene, topology_fluoroethylene):
    maker = TopologyMaker(geometry_fluoroethylene)

    # check every bond was found
    for bond in topology_fluoroethylene.bonds:
        assert bond in maker.g.edges
