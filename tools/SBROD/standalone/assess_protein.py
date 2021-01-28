# coding: utf-8

import os
import sys
import glob
import argparse
import numpy as np
import sys
from sys import exit
from protein_assessor import ProteinAssessor


__author__ = 'Mikhail Karasikov and Sergei Grudinin'
__date__ = '10/03/2017'
__version__ = 1.0
__description__ = 'Smooth Backbone-Reliant Orientation-Dependent Protein Quality Assessment'


BASE_DIR = os.path.dirname(sys.executable) + '/base/'
MODEL = 'pipeline.pkl'
EXECUTABLE = 'Rotamers'
DEFAULT_SIGMA = 0.1866


def restricted_float(x):
    x = float(x)
    if x < 0.0:
        raise argparse.ArgumentTypeError("%r is negative"%(x,))
    return x


def main(pdb_patterns, sigma, transform=False):
    try:
        assessor = ProteinAssessor(BASE_DIR + MODEL, BASE_DIR + EXECUTABLE, sigma)
    except FileNotFoundError as e:
        sys.stderr.write("Error: ranking model '%s' was not found in %s\nAdd ranking model to %s\n" % (MODEL, BASE_DIR, BASE_DIR))
        exit(1)
    except Exception as e:
        sys.stderr.write(str(e) + '\n')
        sys.stderr.write("Error: problem with feature generator. Reinstall the program\n")
        exit(1)

    for pdb_pattern in pdb_patterns:
        pdb_files = glob.glob(pdb_pattern)

        for pdb_file in pdb_files:
            score, err = assessor.predict_from_pdb(pdb_file)
            sys.stderr.write(err)
            print('{}\t{}'.format(
                pdb_file,
                1 / (1 + np.exp(-5 * (score - 1.6))) if transform and score is not None else score
            ))


if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser(
            description=__description__ + '\nDeveloped by ' + __author__ + ', ' + __date__,
            epilog='This program performs ranking of the protein models.\n'
                   'The larger the output score for a protein model, the better quality for that protein model is predicted.\n'
                   'By default, the predicted scores fall into (-infty, +infty).\n'
                   'Contact the authors: karasikov@phystech.edu',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'pdb',
            nargs='+',
            help='pdb files with protein structures (eg ../dataset/*.pdb)'
        )
        parser.add_argument(
            '--sigma',
            type=restricted_float,
            default=DEFAULT_SIGMA,
            help='smoothing parameter for feature generation'
        )
        parser.add_argument(
            '--scale',
            action='store_true',
            help='transform scores to interval (0,1) preserving the ranking'
        )
        parser.add_argument(
            '-v', '--version',
            action='version',
            version='%(prog)s ' + str(__version__)
        )
        args = parser.parse_args()

        main(args.pdb, args.sigma, args.scale)
    except KeyboardInterrupt:
        exit(1)
    except Exception as e:
        sys.stderr.write('Unknown error' + (': ' if str(e) else '') + str(e) + '\n')
        exit(1)
