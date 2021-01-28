#sbrod
import os
import glob
import numpy as np
import scipy
from scipy.io import loadmat
import pandas as pd
import pickle
import hashlib


def dataset_binary_filename(mat_file_pattern, checksum, binaries_path):
    return binaries_path + '/' + mat_file_pattern + checksum


def aggregate_dataset(path_to_dataset, mat_file_pattern='*.mat', binaries_path='binaries'):
    """Loads sparse features and labels"""

    X, scores_records = [], []

    protein_dirs = glob.glob(path_to_dataset + '/*')
    print('Subfolders in the directory: {}'.format(len(protein_dirs)))

    loaded_features_list = []
    for protein_dir in sorted(protein_dirs):
        structure_paths = glob.glob(protein_dir + '/' + mat_file_pattern)
        if len(structure_paths) == 0:
            continue

        scores = pd.read_table(protein_dir + '/scores.txt', index_col='NAME')
        scores.index = protein_dir + '/' + scores.index

        for structure_path in structure_paths:
            try:
                x = loadmat(structure_path)
            except:
                print('\nParsing features error:', structure_path)
                continue
            try:
                x = x['SpMat_' + structure_path.split('/')[-1].split('.')[-2]].T
                structure_name = '*'.join(structure_path.split('/')[-1].split('*')[:-1])
                scores_records.append(
                    scores.loc[protein_dir + '/' + structure_name].astype(float)
                )
                X.append(x)
                loaded_features_list.append(structure_path)
            except:
                print('\nName or scores parsing error:', structure_path)
                continue

    if len(loaded_features_list) == 0:
        raise Exception('Empty dataset')

    X, scores = scipy.sparse.vstack(X), pd.DataFrame(scores_records)

    checksum_hash = hashlib.md5(''.join(scores.index.sort_values()).encode()).hexdigest()
    print('md5 hash: {}'.format(checksum_hash))

    pickle_dataset(dataset_binary_filename(mat_file_pattern, checksum_hash, binaries_path), X, scores)
    print('Dataset was pickled')

    for loaded_features in loaded_features_list:
        os.remove(loaded_features)


def save_sparse_csr(filename, array):
    np.savez_compressed(filename, data=array.data, indices=array.indices,
                        indptr=array.indptr, shape=array.shape)


def pickle_dataset(name, X, y):
    save_sparse_csr(name, X)
    y.to_csv(name + '.csv')
