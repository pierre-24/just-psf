"""
"Just get me a topology, for god’s sake!"
Create a topology, based on the distance matrix.
"""

import argparse
import sys

from just_psf.geometry import Geometry
from just_psf.typology_maker import TopologyMaker


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', type=argparse.FileType('r'), help='input geometry')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='output', default=sys.stdout)

    args = parser.parse_args()

    # read files
    geometry = Geometry.from_xyz(args.infile)

    # make topology
    TopologyMaker(geometry).topology().to_psf(args.output)


if __name__ == '__main__':
    main()
