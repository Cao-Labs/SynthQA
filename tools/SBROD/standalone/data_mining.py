import numpy as np
from scipy import sparse as sp
from scipy.io import loadmat
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression, Ridge


"""Data loading"""


def csr_vstack(a, b):
    """Takes 2 matrices and appends the second one to the bottom of the first one"""
    assert(type(a) is sp.csr_matrix)
    assert(type(b) is sp.csr_matrix)
    a.data = np.hstack((a.data, b.data))
    a.indices = np.hstack((a.indices, b.indices))
    a.indptr = np.hstack((a.indptr, (b.indptr.astype(np.int64) + a.nnz)[1:]))
    a._shape = (a.shape[0] + b.shape[0], b.shape[1])
    return a


class OneHotImputer(BaseEstimator, TransformerMixin):
    def __init__(self, weight=1):
        self.weight = weight

    def fit(self, X, scores):
        raise Exception('Use fit_transform')
        return self

    def fit_transform(self, X, scores):
        one_hot_labels = self.weight * one_hot_from_scores(scores)
        self.n = one_hot_labels.shape[1]
        return csr_vstack(X.T.tocsr(), sp.csc_matrix(one_hot_labels).T).T.tocsr()

    def transform(self, X, scores=None):
        return csr_vstack(X.T.tocsr(), sp.csc_matrix((X.shape[0], self.n)).T).T.tocsr()


class CombinedScaler(BaseEstimator, TransformerMixin):
    def __init__(self, normalizers, nonzero_features=None):
        self.normalizers = normalizers
        self.nonzero_features = nonzero_features

    def fit(self, X, scores):
        raise Exception('Not allowed')

    def fit_transform(self, X, scores):
        raise Exception('Not allowed')

    def transform(self, X, scores=None):
        X = sp.hstack(
            [normalizer.transform(x) for (x, normalizer) in zip(X, self.normalizers)]
        ).tocsc()
        if self.nonzero_features:
            X = X[:, self.nonzero_features]
        return X.tocsr()


class RRModel(Ridge):
    def fit(self, X, scores, **fit_params):
        super(RRModel, self).fit(X, 1 - scores['GDT-TS-score'], **fit_params)
        return self
