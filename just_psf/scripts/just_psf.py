"""
"Just get me a topology, for god’s sake!"
Create a Protein Structure File (PSF), based on the distance matrix.
"""

import argparse
import sys

from just_psf.geometry import Geometry
from just_psf.geometry_analyzer import GeometryAnalyzer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', type=argparse.FileType('r'), help='input geometry')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='output', default=sys.stdout)

    args = parser.parse_args()

    # read file
    geometry = Geometry.from_xyz(args.infile)

    # make topology
    # "ext xplor" format required, because atom types may be longer than 4 chars!
    GeometryAnalyzer(geometry).structure().to_psf(args.output, flags=['EXT', 'XPLOR'])


if __name__ == '__main__':
    main()
