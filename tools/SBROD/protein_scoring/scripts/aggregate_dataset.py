import os
import glob
import numpy as np
import scipy
from scipy.io import loadmat
#import pandas as pd
import pickle
import hashlib

from sklearn.decomposition import TruncatedSVD


def dataset_binary_filename(mat_file_pattern, checksum, binaries_path):
    return binaries_path + '/' + mat_file_pattern + checksum


def aggregate_dataset(path_to_dataset, mat_file_pattern='*.mat', binaries_path='binaries', working_folder='./'):
    """Loads sparse features and labels"""

    X, scores_records = [], []

    protein_dirs = [working_folder + "/mat_files/"+path_to_dataset.split('/')[-1]] #glob.glob(path_to_dataset + '/*')
    print("Protein directory: %s" % protein_dirs)
    #print('Subfolders in the directory: {}'.format(len(protein_dirs)))

    loaded_features_list = []
    loaded_structure_names = []
    for protein_dir in sorted(protein_dirs):
        print("\n################ ", mat_file_pattern, " ################")

        structure_paths = glob.glob(protein_dir + '/' + mat_file_pattern)
        if len(structure_paths) == 0:
            continue

        #scores = pd.read_table(protein_dir + '/scores.txt', index_col='NAME')
        #scores.index = protein_dir + '/' + scores.index

        for structure_path in structure_paths:
            try:
                x = loadmat(structure_path)
            except:
                print('\nParsing features error:', structure_path)
                continue
            try:
                x = x['SpMat_' + structure_path.split('/')[-1].split('.')[-2]].T
                structure_name = structure_path.split('/')[-1].rstrip(mat_file_pattern[1:]).split('.')[0] #'*'.join(structure_path.split('/')[-1].split('*')[:-1])
                #scores_records.append(
                #    scores.loc[protein_dir + '/' + structure_name].astype(float)
                #)
                X.append(x)
                loaded_features_list.append(structure_path)
                loaded_structure_names.append(structure_name)
                print("Successfully Processed: ", structure_path)
            except:
                print('\nName or scores parsing error:', structure_path)
                continue

    if len(loaded_features_list) == 0:
        raise Exception('Empty dataset')

    # If there's only one input model, duplicate its data so that it works with SVD
    if len(X) == 1:
        X.append(X[0])

    X = scipy.sparse.vstack(X)

    # Perform decomposition
    svd = TruncatedSVD(n_components=min(X.shape)-1, algorithm="arpack") # may have problems if there's only one pdb
    X_new = svd.fit_transform(X)
    #print(X_new, "\n", X_new.shape)

    # Pick the column with the most explained variance
    variances = svd.explained_variance_ratio_
    print("--- Explained Variance: ", variances)
    col = 0
    if len(variances) > 1: # if there's only one column, the variance can be nan or inf
        col = variances.tolist().index( max(variances) )

    # Pair the values with the file name
    final_feature_data = []
    for i in range(len(loaded_structure_names)):
        model_filename = loaded_structure_names[i]
        value = X_new[i][col]
        final_feature_data.append( (model_filename, value) )

    print('\n\n')

    subfolder_name = protein_dirs[0].split('/')[-1]
    final_out_dir = binaries_path + '/' + subfolder_name
    if not os.path.isdir(final_out_dir):
        os.makedirs(final_out_dir)

    output_path =  final_out_dir + '/' + mat_file_pattern[2:-4] + ".sbrod"
    with open(output_path, 'wt') as outfile:
        for entry in final_feature_data:
            row = entry[0] + " " + str(entry[1]) + "\n"
            print(row.rstrip('\n'))
            outfile.write(row)

    #checksum_hash = hashlib.md5(''.join(scores.index.sort_values()).encode()).hexdigest()
    #print('md5 hash: {}'.format(checksum_hash))

    #pickle_dataset(dataset_binary_filename(mat_file_pattern, checksum_hash, binaries_path), X, scores)
    #print('Dataset was pickled')

    #for loaded_features in loaded_features_list:
       # os.remove(loaded_features)


def save_sparse_csr(filename, array):
    np.savez_compressed(filename, data=array.data, indices=array.indices,
                        indptr=array.indptr, shape=array.shape)


def pickle_dataset(name, X, y):
    save_sparse_csr(name, X)
    y.to_csv(name + '.csv')
