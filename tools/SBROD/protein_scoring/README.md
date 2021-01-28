# Protein Quality Assessment

Source code for training <b>S</b>mooth <b>B</b>ackbone-<b>R</b>eliant <b>O</b>rientation-<b>D</b>ependent protein quality assessment.

## Typical workflow for training a protein quality assessment method
1. Prepare data (see [datasets](../datasets))
2. [Extract features](#feature-extraction) (or download [features already extracted](#download-features))
3. Train a scoring model
4. Evaluate, compare to other methods

### Feature extraction

`SBROD` uses geometric features, which are extracted by the tool `Rotamers`.

`Rotamers` is written in C++ for high efficiency and used for modeling spacial protein structures and various other tasks.

Its standalone binary executables are included in the `SBROD` package.
1. Download and unpack standalone executables from [here](../README.md#download-precompiled-standalone-executables)
2. Extract and copy `Rotamers` from `sbrod/base/` to [`scripts/`](./scripts/)
3. Run [`scripts/generate_dataset.sh`](./scripts/generate_dataset.sh) and see other [scripts](./scripts)

Scripts in [scripts](./scripts) assume that the data is organized in the following directory structure:
* dataset_1
* ... (datasets, e.g. CASP11, CASP12)
* dataset_N
    * target_dir_1
    * ... (each target_dir contains predictions for one native structure, or primary sequence)
    * target_dir_K
        - protein_model_1
        - ... (predicted protein structures, may include the native structure)
        - protein_model_M
        - `scores.txt` (file with computed ground truth, see [compute-ground-truth-for-decoy-protein-models](../datasets/README.md#compute-ground-truth-for-decoy-protein-models))

#### Download features
To download preprocessed data with features already extracted from various `CASP` datasets, use [this link](https://drive.google.com/open?id=1uh8Zv0-ZDZLLxAYQu5dhwi-icQ5mRflZ).

#### Feature parameters
All developed geometric features are parametric functions of protein structures.
To tune the parameters and change feature extraction procedure see [gridsearch](./gridsearch).

### Training a scoring model
Once features are extracted, one can train scoring models.

See various scoring models proposed in [training.ipynb](./training.ipynb)
or train a default linear model ([train_final_model.ipynb](./train_final_model.ipynb) or [train_all.py](./train_all.py)).


### Evaluation
`SBROD` was evaluated against other methods for single-model protein quality assessment (ProQ2, ProQ3, VoroMQA, RWplus).

Find the results in [comparison/CASP12.ipynb](./comparison/CASP12.ipynb) and [other notebooks](./comparison).

Follow the instructions in [corresponding folders](./comparison/others) to reproduce those results or perform a new comparison on a different dataset.

## Experiments

Find other experiments in [experiments](./experiments).
