# coding: utf-8

import os
import re
import joblib
import subprocess
import numpy as np
from scipy.io import loadmat
from tempfile import NamedTemporaryFile
import data_mining


__author__ = 'Mikhail Karasikov'


class ProteinAssessor():
    def __init__(self, path_to_scorer, executable, sigma):
        self.executable = executable
        try:
            subprocess.check_output([self.executable, '--help'], stderr=subprocess.STDOUT)
        except Exception as e:
            raise Exception(e)

        self.params_list = [
            '--featurization residues -d 4 -b 10 -a 12 -c 5 -n 0 -s {} --skip_errors'.format(sigma),
            '--featurization hbonds -b 6 -a 6 -c 6 -n 2 -s {} --skip_errors'.format(sigma),
            '--featurization solvation -b 3 -a 2 -c 15 -s {} --skip_errors'.format(sigma),
            '--featurization backboneatom -b 25 -c 7 -n 0 -s {} --residue_type_dependent --skip_errors'.format(sigma)
        ]
        self.scorer = joblib.load(path_to_scorer)

    def predict_from_pdb(self, pdb_filename):
        try:
            with open(pdb_filename, 'r') as f_in:
                return self.__predict(f_in)    
        except Exception as e:
            return None, "Error"

    def predict_from_str(self, string):
        try:
            f_in = NamedTemporaryFile('w', bufsize=100, suffix='.pdb')
            f_in.write(string)
            return self.__predict(f_in)
        except Exception as e:
            return None, "Error"

    def __predict(self, f_in):
        mat_filenames, output = self.__featurize(f_in)
        if mat_filenames is None:
            return None, output

        features = self.__load_features(mat_filenames)
        return self.__score(features), output

    def __load_features(self, mat_filenames):
        features = []
        for mat_filename in mat_filenames:
            x = loadmat(mat_filename)
            os.remove(mat_filename)
            x = x['SpMat_' + re.split('[(\\\\)|/]', mat_filename)[-1].split('.')[-2]]
            features.append(x.T)
        return features

    def __featurize(self, pdb_file):
        mat_filenames = []
        output = []
        for params in self.params_list:                
            mat_file = NamedTemporaryFile(suffix='.mat')
            mat_file.close()

            try:
                subprocess.check_output(
                    [self.executable, '--mode', 'featurize'] + params.split() + ['-i', pdb_file.name, '-o', mat_file.name[:-4]],
                    stderr=subprocess.STDOUT
                ).decode()
                mat_filenames.append(mat_file.name)
            except subprocess.CalledProcessError as e:
                output.append(e.output)
                for mat_filename in mat_filenames:
                    os.remove(mat_filename)
                return None, '\n'.join(output)

        return mat_filenames, '\n'.join(output)

    def __score(self, x):
        return 1 - float(self.scorer.predict(x))

    def complete(self, string):
        try:
            pdb_in = NamedTemporaryFile('w', bufsize=100, suffix='.pdb')
            pdb_in.write(string)

            pdb_out = NamedTemporaryFile(suffix='.pdb')

            try:
                subprocess.check_output(
                    [self.executable, '--mode', 'complete', '-i', pdb_in.name, '-o', pdb_out.name],
                    stderr=subprocess.STDOUT
                ).decode()
            except subprocess.CalledProcessError as e:
                return None, e.output

            completed_pdb = pdb_out.read()

            return completed_pdb, ''

        except Exception as e:
            return None, "Error"
