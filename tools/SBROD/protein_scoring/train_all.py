import numpy as np
import scipy
import pandas as pd
import joblib
from sklearn.preprocessing import Normalizer, MaxAbsScaler

from data_mining import *


residues_normalizer = Normalizer()
hbonds_normalizer = MaxAbsScaler()
solvation_normalizer = MaxAbsScaler()
backboneatom_normalizer = Normalizer()

CASP = [
    ('*residues-d4-b10-a12-c5-n0--skip_errors.mat',
     '961920c35dbee45086ef2544978c884f',
     lambda X: residues_normalizer.fit_transform(X),
     lambda X: residues_normalizer.transform(X)),
    ('*hbonds-b6-a6-c6-n2--skip_errors.mat',
     '961920c35dbee45086ef2544978c884f',
     lambda X: hbonds_normalizer.fit_transform(X),
     lambda X: hbonds_normalizer.transform(X)),
    ('*solvation-b3-a2-c15--skip_errors.mat',
     '961920c35dbee45086ef2544978c884f',
     lambda X: solvation_normalizer.fit_transform(X),
     lambda X: solvation_normalizer.transform(X)),
    ('*backboneatom-b25-c7-n0--residue_type_dependent--skip_errors.mat',
     '6ffffc5025bb8a4c3ad76dcd0d823b1e',
     lambda X: backboneatom_normalizer.fit_transform(X),
     lambda X: backboneatom_normalizer.transform(X)),
]

NMA = [
    ('*residues-d4-b10-a12-c5-n0--skip_errors.mat',
     'a34560ffd5d69bdce838fbe729b84525',
     lambda X: residues_normalizer.fit_transform(X),
     lambda X: residues_normalizer.transform(X)),
    ('*hbonds-b6-a6-c6-n2--skip_errors.mat',
     'a34560ffd5d69bdce838fbe729b84525',
     lambda X: hbonds_normalizer.fit_transform(X),
     lambda X: hbonds_normalizer.transform(X)),
    ('*solvation-b3-a2-c15--skip_errors.mat',
     'a34560ffd5d69bdce838fbe729b84525',
     lambda X: solvation_normalizer.fit_transform(X),
     lambda X: solvation_normalizer.transform(X)),
    ('*backboneatom-b25-c7-n0--residue_type_dependent--skip_errors.mat',
     'a34560ffd5d69bdce838fbe729b84525',
     lambda X: backboneatom_normalizer.fit_transform(X),
     lambda X: backboneatom_normalizer.transform(X)),
]

CASP_predictions = [
    ('*residues-d4-b10-a12-c5-n0--skip_errors.mat',
     'c7d3ecb2312adf488b7fc9081c0f6bf5',
     lambda X: residues_normalizer.fit_transform(X),
     lambda X: residues_normalizer.transform(X)),
    ('*hbonds-b6-a6-c6-n2--skip_errors.mat',
     'c7d3ecb2312adf488b7fc9081c0f6bf5',
     lambda X: hbonds_normalizer.fit_transform(X),
     lambda X: hbonds_normalizer.transform(X)),
    ('*solvation-b3-a2-c15--skip_errors.mat',
     'c7d3ecb2312adf488b7fc9081c0f6bf5',
     lambda X: solvation_normalizer.fit_transform(X),
     lambda X: solvation_normalizer.transform(X)),
    ('*backboneatom-b25-c7-n0--residue_type_dependent--skip_errors.mat',
     'c7d3ecb2312adf488b7fc9081c0f6bf5',
     lambda X: backboneatom_normalizer.fit_transform(X),
     lambda X: backboneatom_normalizer.transform(X)),
]


X_CASP, scores_CASP = get_dataset(
    [(pattern, checksum, fit_transform)
         for (pattern, checksum, fit_transform, transform) in CASP],
    '^.*CASP([5-9]|10|11|12)/T..../.*$'
)
scores_CASP = scores_CASP[['RMSD', 'GDT-TS-score']]

X_NMA, scores_NMA = get_dataset(
    [(pattern, checksum, transform)
         for (pattern, checksum, fit_transform, transform) in NMA],
    '^.*CASP([5-9]|10|11|12)/T..../.*$'
)
scores_NMA = scores_NMA[['RMSD', 'GDT-TS-score']]

X_CASP_predictions, scores_CASP_predictions = get_dataset(
    [(pattern, checksum, transform)
         for (pattern, checksum, fit_transform, transform) in CASP_predictions],
    '^.*CASP([5-9]|10|11|12)Predictions/T..../.*$'
)
scores_CASP_predictions.rename(columns={'GDT_TS': 'GDT-TS-score'}, inplace=True)
scores_CASP_predictions['GDT-TS-score'] /= 100
scores_CASP_predictions = scores_CASP_predictions[['RMSD', 'GDT-TS-score']]

X, scores = combine_datasets(
    (X_CASP, scores_CASP),
    (X_NMA, scores_NMA),
    (X_CASP_predictions, scores_CASP_predictions)
)

del X_CASP, scores_CASP, X_NMA, scores_NMA, X_CASP_predictions, scores_CASP_predictions

def dump_model(model, filename):
    joblib.dump(
        Pipeline([('scaler', CombinedScaler([residues_normalizer, hbonds_normalizer,
                                             solvation_normalizer, backboneatom_normalizer])),
                  ('scorer', model)]),
        filename, protocol=2
    )


ridge_pipeline = Pipeline([
    ('imputer', OneHotImputer(100)),
    ('model', RRModel(copy_X=False, normalize=False, fit_intercept=False, solver='sparse_cg', alpha=5))
])

ridge_pipeline.fit(X, scores)

dump_model(ridge_pipeline, 'ridge_pipeline_{}_{}.pkl'.format(*X.shape))

print("Done")
