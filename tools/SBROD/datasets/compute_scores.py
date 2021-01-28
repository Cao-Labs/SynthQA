#!/usr/bin/python2.7

""""
This script measures the structural similarity between predicted
protein models and their corresponding native structures.
(A wrapper for running TMscore in multiple threads)

Dependency: the TMscore tool https://zhanglab.ccmb.med.umich.edu/TM-score

Usage: ./compute_scores.py <structures_paths> <native_pattern> <decoy_pattern> <num_threads>

Example: ./compute_scores.py "CASP/data/CASP10/*" "T0*.pdb" "*" 8
"""

import sys
from glob import glob
import os
import commands
from multiprocessing import Pool


scores_filename = "scores.txt"


def on_error(decoy_file):
    os.remove(decoy_file)
    print("%s was removed" % decoy_file)


def compute_scores_for_pair(first_filename, second_filename):
    output = commands.getstatusoutput("./TMscore {} {}".format(first_filename, second_filename))
    if output[1].startswith(' There is no common residues in the input structures'):
        output = commands.getstatusoutput("./TMscore -c {} {}".format(first_filename, second_filename))
        if output[1].startswith(' There is no common residues in the input structures'):
            raise Exception(output[1] + "\nReturn value: {}".format(output[0]))

    if output[0] or output[1].startswith('Warning'):
        raise Exception(output[1] + "\nReturn value: {}".format(output[0]))

    try:
        lines = output[1].split('\n')
        rmsd = float(lines[14][29:])
        tm = float(lines[16][13:21])
        maxsub = float(lines[17][13:21])
        gdtts = float(lines[18][13:21])
        gdtha = float(lines[19][13:21])
        return {
            "RMSD": rmsd,
            "TM-score": tm,
            "MaxSub-score": maxsub,
            "GDT-TS-score": gdtts,
            "GDT-HA-score": gdtha
        }

    except:
        raise Exception(output[1] + "\nReturn value: {}".format(output[0]))


def compute_scores_for_decoy_set(args):
    target, native_pattern, decoy_pattern = args

    print target
    if os.path.isfile(target):
        return

    native_file = target + "/" + native_pattern
    decoys_glob = target + "/" + decoy_pattern

    with open(target + "/" + scores_filename, "w") as listfile:
        listfile.write(
            "NAME" +
            "\tRMSD\tTM-score\tMaxSub-score\tGDT-TS-score\tGDT-HA-score" +
            "\tRMSD-backwards\tTM-score-backwards\tMaxSub-score-backwards\tGDT-TS-score-backwards\tGDT-HA-score-backwards\n"
        )

        # loop through the mutants
        for decoy_file in glob(decoys_glob):
            if decoy_file.endswith(scores_filename):
                continue

            try:
                scores = compute_scores_for_pair(decoy_file, native_file)
                scores_backwards = compute_scores_for_pair(native_file, decoy_file)
                decoy_name = decoy_file.split("/")[-1]
            except Exception as e:
                print("Exception when processing structure: %s" % decoy_file)
                print(e)
                on_error(decoy_file)
                continue

            listfile.write(
                decoy_name +
                "\t{RMSD}\t{TM-score}\t{MaxSub-score}\t{GDT-TS-score}\t{GDT-HA-score}".format(**scores) +
                "\t{RMSD}\t{TM-score}\t{MaxSub-score}\t{GDT-TS-score}\t{GDT-HA-score}\n".format(**scores_backwards)
            )


def main():
    if len(sys.argv) != 5:
        print("Usage: {} <structures_paths> <native_pattern> <decoy_pattern> <num_threads>".format(sys.argv[0]))
        exit(1)

    _, structures_paths, native_pattern, decoy_pattern, num_threads = sys.argv

    targets = glob(structures_paths)

    pool = Pool(processes=int(num_threads))
    pool.imap_unordered(
        compute_scores_for_decoy_set,
        [(target, native_pattern, decoy_pattern) for target in targets]
    )
    pool.close()
    pool.join()

    print "finished!"


if __name__ == '__main__':
    main()
