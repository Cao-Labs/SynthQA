#!/usr/bin/python3
#sbrod
"""
Parallel feature extraction

Dependency: Rotamers (the SBROD feature extractor)
https://team.inria.fr/nano-d/software/sbrod/
Download the SBROD method and extract Rotamers from the `base` directory

Example: python featurize_dataset_parallel.py "../../datasets/CASP/data/CASP*/*/*[!t]" backboneatoms "-b 10" 8
"""

import glob
import sys
import subprocess
import re
from multiprocessing.dummy import Pool


ROTAMERS_EXECUTABLE = "./Rotamers"


def generate_features(args):
    featurization_type, file, args_str = args

    print(file)
    completed_process = subprocess.run([ROTAMERS_EXECUTABLE, "--mode", "featurize",
                                                             "--featurization", featurization_type,
                                        "-i", file,
                                        "-o", file + '*' + featurization_type + args_str.replace(' ', '')] +
                                        args_str.split(' '))
    if completed_process.returncode:
        print("ERROR occured")
        #print(completed_process.stderr)
        #rm -r $(dirname $file)


def main():
    if len(sys.argv) != 5:
        print('Usage: {} <pdb_pattern> <featurization_type> <arguments> <num_processes>'.format(sys.argv[0]))
        exit(1)


    pdb_pattern, featurization_type, args_str, num_processes = sys.argv[1:]

    pdb_files = glob.glob(pdb_pattern)

    pool = Pool(processes=int(num_processes))
    pool.map(generate_features, ((featurization_type, file, args_str) for file in pdb_files))


if __name__ == '__main__':
    main()
