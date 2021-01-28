import os
import glob
import numpy as np
import scipy
from scipy.io import loadmat
import pandas as pd
from scipy.stats import pearsonr, spearmanr, kendalltau
from matplotlib import pyplot as plt
from sklearn.feature_extraction import DictVectorizer
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold
import pickle
from itertools import combinations
from IPython.display import display
import re
import itertools
import hashlib
import seaborn


"""Data loading"""

def reglob(path, exp, invert=False):
    m = re.compile(exp)

    if invert:
        return [f for f in glob.glob(path) if not m.match(f)]

    return [f for f in glob.glob(path) if m.match(f)]


def dataset_binary_filename(mat_file_pattern, checksum, binaries_path):
    return binaries_path + '/' + mat_file_pattern + checksum


def load_dataset(path_to_dataset, num_proteins=np.inf, mat_file_pattern='*.d4_b10.mat',
                 decoy_dirs_pattern='.*', checksum='', binaries_path='binaries'):
    """Loads sparse features and labels"""

    try:
        if checksum is None:
            raise Exception()
        X, scores = load_pickled_dataset(dataset_binary_filename(mat_file_pattern, checksum, binaries_path))
    except:
        X, scores_records = [], []

        features_to_load = {}
        protein_dirs = glob.glob(path_to_dataset + '/*')
        print('Subfolders in the directory: {}'.format(len(protein_dirs)))
        for protein_dir in protein_dirs[:min(len(protein_dirs), num_proteins)]:
            features_to_load[protein_dir] = glob.glob(protein_dir + '/' + mat_file_pattern)

        loaded_features_list = []
        for protein_dir, structure_paths in sorted(features_to_load.items()):
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
                    x = x['SpMat_' + structure_path.split('/')[-1].split('.')[-2]]
                    feature_vector = x.T
                    structure_name = '*'.join(structure_path.split('/')[-1].split('*')[:-1])
                    scores_records.append(
                        scores.loc[protein_dir + '/' + structure_name].astype(float)
                    )
                    X.append(feature_vector)
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

        if checksum_hash == checksum or checksum is None:
            for loaded_features in loaded_features_list:
                os.remove(loaded_features)
        else:
            print('Different checksums')

    return select_subset(X, scores, decoy_dirs_pattern)


def select_subset(X, scores, decoy_dirs_pattern):
    m = re.compile(decoy_dirs_pattern)
    indices = np.array([bool(m.match(f)) for f in scores.index])
    return X[indices], scores[indices].copy()


def one_hot_from_scores(scores):
    names = ['/'.join(structure_path.split('/')[:-1]) for structure_path in scores.index]
    return DictVectorizer(sparse=False).fit_transform([{name: 1} for name in names])


def save_sparse_csr(filename, array):
    np.savez_compressed(filename, data=array.data, indices=array.indices,
                        indptr=array.indptr, shape=array.shape)


def load_sparse_csr(filename):
    loader = np.load(filename)
    return scipy.sparse.csr_matrix(
        (loader['data'], loader['indices'], loader['indptr']),
        shape=loader['shape']
    )


def load_pickled_dataset(name):
    X = load_sparse_csr(name + '.npz')
    y = pd.read_csv(name + '.csv', header=0, index_col=0)
    return X, y


def pickle_dataset(name, X, y):
    save_sparse_csr(name, X)
    y.to_csv(name + '.csv')


def remove_zero_features(X):
    s = X.sum(0)
    nonzero_features = s.nonzero()[1]
    return X[:, nonzero_features].tocsr(), nonzero_features


def generate_benchmark(results, performance_quality):
    benchmark = pd.DataFrame()
    for mat_file_pattern in results.keys():
        for model_name in results[mat_file_pattern].keys():
            benchmark.ix[mat_file_pattern, model_name] = performance_quality(results[mat_file_pattern][model_name])

    for i in range(benchmark.shape[0]):
        params = '.'.join(benchmark.index[i].split('.')[:-1])

        for parameter, value in zip([k[1:] for k in re.findall('-[a-z_]+', params)],
                                    re.split('[-]+[a-z_]+', params)[1:]):
            benchmark.ix[i, parameter] = value if value != '' else 1

    benchmark = benchmark.astype(float)

    return benchmark


"""Evaluation and Visualization"""

def plot_heatmaps(benchmark, values, parameters, info='', num_cols=3, figsize=0.8):
    benchmark = benchmark.copy()

    for index, columns in combinations(parameters, 2):
        if len(parameters) == 2:
            benchmark[''] = ''
            parameters = list(parameters) + ['']

        groups = list(set(parameters) - set([index, columns]))
        groups_values = sorted(set([tuple(r) for r in benchmark[groups].values]))

        num_rows = int(np.ceil(len(groups_values) / num_cols))

        fig, ax = plt.subplots(num_rows, num_cols, figsize=(num_cols * 4 * figsize, num_rows * 3.5 * figsize))
        for i, groups_value in enumerate(groups_values):
            benchmark_slice = benchmark[(benchmark[groups].values == groups_value).all(1)]
            pivot_table = pd.pivot_table(benchmark_slice, values=values, columns=[columns], index=[index])
            pivot_table.index.name = '${}$'.format(pivot_table.index.name)
            pivot_table.columns.name = '${}$'.format(pivot_table.columns.name)

            seaborn.heatmap(pivot_table, ax=ax.item(i), vmin=benchmark[values].min(),
                                                        vmax=benchmark[values].max())

            parameters_title = ', '.join(['{}={}'.format(group, value) for (group, value) in zip(groups, groups_value)])
            ax.item(i).set_title('${}$'.format(parameters_title), size=15)
            ax.item(i).xaxis.label.set_size(15)
            ax.item(i).yaxis.label.set_size(15)

        for i in range(i + 1, num_rows * num_cols):
            ax.item(i).axis('off')

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.suptitle(', '.join([values, info]), size=16)
        plt.show()


def plot_correlations(pairs, block_names, num_cols=3, score_name='Actual score', find_native=np.argmin):
    num_rows = int(np.ceil(len(pairs) / num_cols))
    fig, ax = plt.subplots(num_rows, num_cols, figsize=(num_cols * 4, num_rows * 4))
    for i, (predicted_score, actual_score) in enumerate(pairs):
        ax.item(i).scatter(actual_score, predicted_score)
        ax.item(i).scatter([actual_score[find_native(actual_score)]],
                           [predicted_score[find_native(actual_score)]], s=50, c='r')
        ax.item(i).set_xlabel(score_name, fontsize=10)
        ax.item(i).set_ylabel('Predicted score', fontsize=10)
        block_corr = safe_correlation(pearsonr, predicted_score, actual_score)[0]
        ax.item(i).set_title('{}, corr={:.3f}'.format(block_names[i], block_corr))
    for i in range(i + 1, num_rows * num_cols):
        ax.item(i).axis('off')
    #fig.suptitle('Correlations', fontsize=16)
    plt.tight_layout()
    plt.show()


def weighted_spearmanr(x, y, f=lambda x: 1 - np.tanh(x / 4)):
    if len(x) != len(y):
        raise Exception('x and y have to be the same length')
    x, y = np.array(x), np.array(y)
    n = len(x)
    x_rank, y_rank = x.argsort().argsort(), y.argsort().argsort()
    weights = f(x)
    return 1 - 2 * sum((x_rank - y_rank) ** 2 * weights /
                       sum((n - 1 - 2 * np.arange(n)) ** 2 * sorted(weights)))


def scoring_results(X, scores, score_function, num_plots=24, plot_score='RMSD'):
    if not (X.shape[0] == len(scores)):
        raise Exception('Dimensions error')

    one_hot_labels = one_hot_from_scores(scores)
    predicted_scores = [score_function(X[block]) for block in one_hot_labels.T.astype(bool)]
    block_names = [scores.index[block][0].split('/')[-2] for block in one_hot_labels.T.astype(bool)]

    results = pd.DataFrame(
        index=['Mean rank of the native',
               'Top 1', 'Top 5',
               'Mean score for top 1', 'Mean loss score',
               'Mean Pearson', 'Mean Spearman', 'Mean Kendall tau',
               'PMCC', 'SMCC',
               'wmPMCC', 'wmSMCC',
               'Mean Z-score'],
        columns=scores.columns
    )

    native_ranks, pairs = {}, {}

    for cur_score in scores.columns:
        native_ranks[cur_score], pairs[cur_score], native_scores = [], [], []
        current_ranks, current_pairs = native_ranks[cur_score], pairs[cur_score]

        for block, predicted_scores_block in zip(one_hot_labels.T.astype(bool), predicted_scores):
            actual_scores_block = np.array(scores.ix[block, cur_score])
            current_pairs.append((predicted_scores_block, actual_scores_block))
            native_idx = actual_scores_block.argmin() if cur_score == 'RMSD' else actual_scores_block.argmax()
            current_ranks.append(predicted_scores_block.argsort().argsort()[native_idx])
            native_scores.append(actual_scores_block[native_idx])

        current_ranks = np.array(current_ranks)
        results.loc['Mean rank of the native', cur_score] = np.mean(current_ranks)
        results.loc['Top 1', cur_score] = np.mean(current_ranks == 0)
        results.loc['Top 5', cur_score] = np.mean(current_ranks < 5)
        top_scores = [actual_scores_block[predicted_scores_block.argmin()]
                        for (predicted_scores_block, actual_scores_block) in current_pairs]
        results.loc['Mean score for top 1', cur_score] = np.mean(top_scores)
        results.loc['Mean loss score', cur_score] = np.mean(np.array(native_scores) - np.array(top_scores))
        results.loc['Mean Pearson', cur_score] = np.mean([safe_correlation(pearsonr, *pair)[0] for pair in current_pairs])
        results.loc['Mean Spearman', cur_score] = np.mean([safe_correlation(spearmanr, *pair)[0] for pair in current_pairs])
        results.loc['Mean Kendall tau', cur_score] = np.mean([safe_correlation(kendalltau, *pair)[0] for pair in current_pairs])

        def fisher_transform(r):
            return np.log((1 + r) / (1 - r)) / 2

        def inverse_fisher_transform(z):
            return (np.exp(z) - np.exp(-z)) / (np.exp(z) + np.exp(-z))

        results.loc['PMCC', cur_score] = safe_correlation(pearsonr, list(itertools.chain(*list(zip(*current_pairs))[0])),
                                                              list(itertools.chain(*list(zip(*current_pairs))[1])))[0]
        results.loc['SMCC', cur_score] = safe_correlation(spearmanr, list(itertools.chain(*list(zip(*current_pairs))[0])),
                                                               list(itertools.chain(*list(zip(*current_pairs))[1])))[0]
        results.loc['wmPMCC', cur_score] = inverse_fisher_transform(
            np.mean([fisher_transform(safe_correlation(pearsonr, *pair)[0]) for pair in current_pairs])
        )
        results.loc['wmSMCC', cur_score] = inverse_fisher_transform(
            np.mean([fisher_transform(safe_correlation(spearmanr, *pair)[0]) for pair in current_pairs])
        )
        results.loc['Mean Z-score', cur_score] = np.mean(
            [(actual_scores_block[predicted_scores_block.argmin()] -
              np.mean(actual_scores_block)) / np.std(actual_scores_block)
             for (predicted_scores_block, actual_scores_block) in current_pairs]
        )


    ranks = np.array(native_ranks[plot_score])

    if num_plots > 0:
        plt.figure(figsize=(7, 4))
        plt.hist(ranks, bins=max(ranks) - min(ranks) + 1)
        #plt.xlim([0, max(max(ranks) - min(ranks) + 1, 10)])
        plt.xlabel('Rank')
        plt.ylabel('Number of native structures')
        plt.title('Correct rate: ${:.1f}\%$'.format(100 * (ranks == 0).sum() / len(ranks)), size=12)
        plt.show()

        plot_correlations(pairs[plot_score][:min(len(pairs[plot_score]), num_plots)],
                          block_names[:min(len(pairs[plot_score]), num_plots)],
                          6, plot_score,
                          find_native=np.argmin if plot_score == 'RMSD' else np.argmax)

    return results, (ranks, pairs)


def safe_correlation(func, first, second):
    return func(first, second)
    std = max(np.std(first), np.std(second))
    return func(first + np.random.randn(len(first)) * std * 1e-10,
                second + np.random.randn(len(second)) * std * 1e-10)

"""Dataset transforms"""

def train_test_split(one_hot_labels, test_ratio=0.3, seed=None):
    n = one_hot_labels.shape[1]
    num_test_blocks = round(test_ratio * n)

    np.random.seed(seed)
    test_blocks = np.random.choice(range(n), num_test_blocks, replace=False)
    train_blocks = np.array(list(set(range(n)) - set(test_blocks)))

    test_idx = one_hot_labels[:, test_blocks].sum(1) > 0
    train_idx = one_hot_labels[:, train_blocks].sum(1) > 0

    return train_blocks, train_idx, test_blocks, test_idx


def domain_transform(X, scores, transform=lambda X, scores: (X, scores)):
    X_domains, y_domains = [], []

    for group in one_hot_from_scores(scores).astype(bool).T:
        X_domain, y_domain = transform(X[group], scores[group])
        X_domains.append(X_domain)
        y_domains.append(y_domain)

    return vstack(X_domains), pd.concat(y_domains)

def mean_domain_transform(X, scores):
    X = scipy.sparse.csr_matrix(X) - scipy.sparse.csr_matrix(X.mean(0))[[0] * X.shape[0]]
    scores = scores - scores.mean()
    return X, scores

def top_domain_transform(X, scores):
    top_idx = np.argmax(get_natives(scores).values)
    X = scipy.sparse.csr_matrix(X) - scipy.sparse.csr_matrix(X[top_idx])[[0] * X.shape[0]]
    scores = scores - scores.iloc[top_idx]
    return X, scores

def pairwise_transform(X, scores):
    try:
        X = X.todense()
    except:
        pass
    X, y = np.array(X), np.array(scores)
    X_diff = X[:, None, :] - X[None, :, :]
    y_diff = y[:, None, :] - y[None, :, :]

    X_pairwise = scipy.sparse.csr_matrix(X_diff.reshape(-1, X.shape[1]))
    scores_pairwise = pd.DataFrame(y_diff.reshape(-1, scores.shape[1]), columns=scores.columns)
    mask = np.array(scores_pairwise['GDT-TS-score'] > 0)
    return X_pairwise[mask], scores_pairwise[mask]


def vstack(matrices):
    assert(type(matrices) == tuple or type(matrices) == list)
    for matrix in matrices:
        assert(type(matrix) is scipy.sparse.csr_matrix)

    data = np.hstack((matrix.data for matrix in matrices))

    indices, indptr = [matrices[0].indices], [matrices[0].indptr]
    for matrix in matrices[1:]:
        indptr.append(matrix.indptr.astype(np.int64)[1:] + sum(len(x) for x in indices))
        indices.append(matrix.indices)

    indices = np.hstack(indices)
    indptr = np.hstack(indptr)

    shape = (sum(matrix.shape[0] for matrix in matrices), matrices[0].shape[1])
    return scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)

def csr_vstack(a, b):
    """Takes 2 matrices and appends the second one to the bottom of the first one"""
    assert(type(a) is scipy.sparse.csr_matrix)
    assert(type(b) is scipy.sparse.csr_matrix)
    a.data = np.hstack((a.data, b.data))
    a.indices = np.hstack((a.indices, b.indices))
    a.indptr = np.hstack((a.indptr, (b.indptr.astype(np.int64) + a.nnz)[1:]))
    a._shape = (a.shape[0] + b.shape[0], b.shape[1])
    return a


def combine_datasets(*datasets):
    X = datasets[0][0]
    for Y, _ in datasets[1:]:
        X = csr_vstack(X, Y)
    scores = pd.concat((arg[1] for arg in datasets), axis=0)
    return X, scores


def concatenate_features(X_first, scores_first, X_second, scores_second):
    if X_first is None and scores_first is None:
        return X_second, scores_second

    print('Shape of the first dataset:  {}'.format(X_first.shape))
    print('Shape of the second dataset: {}'.format(X_second.shape))

    def get_indices(scores_first):
        df = scores_first[[]].copy()
        df['i'] = range(df.shape[0])
        return df

    indices = pd.merge(get_indices(scores_first),
                       get_indices(scores_second),
                       left_index=True, right_index=True)

    indices_first = np.array(indices['i_x'])
    indices_second = np.array(indices['i_y'])

    if not np.isclose(scores_first.iloc[indices_first],
                      scores_second.iloc[indices_second], rtol=1e-3).all():
        raise Exception("Scores are not consistent")

    X_first = X_first[indices_first].tocsc().T
    X_second = X_second[indices_second].tocsc().T

    X_first = csr_vstack(X_first, X_second).T
    scores_first = scores_first.iloc[indices_first]
    del X_second, scores_second

    print('Shape of the final dataset:  {}'.format(X_first.shape))

    return X_first, scores_first


class OneHotImputer(BaseEstimator, TransformerMixin):
    def __init__(self, weight=1):
        self.weight = weight

    def fit(self, X, scores):
        raise Exception('Use fit_transform')
        return self

    def fit_transform(self, X, scores):
        one_hot_labels = self.weight * one_hot_from_scores(scores)
        self.n = one_hot_labels.shape[1]
        return csr_vstack(X.T.tocsr(), scipy.sparse.csc_matrix(one_hot_labels).T).T.tocsr()

    def transform(self, X, scores=None):
        return csr_vstack(X.T.tocsr(), scipy.sparse.csc_matrix((X.shape[0], self.n)).T).T.tocsr()


def get_dataset(features, decoy_dirs_pattern):
    def load_features(mat_file_pattern, checksum, decoy_dirs_pattern):
        X, scores = load_pickled_dataset(dataset_binary_filename(mat_file_pattern, checksum, 'binaries'))
        return select_subset(X, scores, decoy_dirs_pattern)

    X, scores = None, None
    for mat_file_pattern, checksum, preprocessing in features:
        X_next, scores_next = load_features(mat_file_pattern, checksum, decoy_dirs_pattern)
        X_next = preprocessing(X_next)
        X, scores = concatenate_features(X, scores, X_next, scores_next)

    return X.tocsr(), scores


class CombinedScaler(BaseEstimator, TransformerMixin):
    def __init__(self, normalizers, nonzero_features=None):
        self.normalizers = normalizers
        self.nonzero_features = nonzero_features

    def fit(self, X, scores):
        raise Exception('Not allowed')

    def fit_transform(self, X, scores):
        raise Exception('Not allowed')

    def transform(self, X, scores=None):
        X = scipy.sparse.hstack(
            [normalizer.transform(x) for (x, normalizer) in zip(X, self.normalizers)]
        ).tocsc()
        if self.nonzero_features:
            X = X[:, self.nonzero_features]
        return X.tocsr()


def get_natives(scores):
    natives = np.array([bool(re.match('(T0...|native)\.pdb', x.split('/')[-1])) for x in scores.index])
    return pd.Series(data=natives, index=scores.index)


class LRModel(LogisticRegression):
    def fit(self, X, scores, **fit_params):
        super(LRModel, self).fit(X, get_natives(scores), **fit_params)
        return self

    def predict(self, X, **predict_params):
        return -super(LRModel, self).decision_function(X)


class RRModel(Ridge):
    def fit(self, X, scores, **fit_params):
        super(RRModel, self).fit(X, 1 - scores['GDT-TS-score'], **fit_params)
        return self


class RRModelNative(RRModel):
    def fit(self, X, scores, **fit_params):
        X, scores = domain_transform(X, scores, top_domain_transform)
        super(RRModelNative, self).fit(X, scores, **fit_params)
        return self


class RankingLR(LogisticRegression):
    def __init__(self, C=1, penalty='l2', fit_intercept=False, rank_transform_threshold=0.1):
        self.rank_transform_threshold=rank_transform_threshold
        super(RankingLR, self).__init__(C=C, penalty=penalty, fit_intercept=fit_intercept)

    def fit(self, X, scores, **fit_params):
        X_rank, scores_rank = domain_transform(X, scores, pairwise_transform)
        mask = np.abs(scores_rank['GDT-TS-score'].values) > self.rank_transform_threshold
        X_rank, scores_rank = X_rank[mask], scores_rank[mask]

        s = (-1) ** np.arange(X_rank.shape[0])
        super(RankingLR, self).fit(np.array(X_rank.todense()) * s.reshape(-1, 1),
                                   scores_rank['GDT-TS-score'].values * s > 0,
                                   **fit_params)
        return self

    def predict(self, X):
        return -super(RankingLR, self).decision_function(X)


def get_folds(scores, n_folds):
    one_hot_labels = one_hot_from_scores(scores)
    n = one_hot_labels.shape[1]

    return ((np.arange(len(scores))[one_hot_labels[:, train].sum(1).astype(bool)],
             np.arange(len(scores))[one_hot_labels[:, test].sum(1).astype(bool)])
                for (train, test) in KFold(n_folds, True, random_state=17).split(np.arange(n))) 


def base_predict(base_model, X):
    return base_model.predict(X)
