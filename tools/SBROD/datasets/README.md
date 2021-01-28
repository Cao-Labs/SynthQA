# Data preparation

This folder contains scripts for downloading and preprocessing data used to train and test protein quality assessment methods.

## Download preprocessed data
All preprocessed datasets with precomputed structure similarity scores are available here:
- [CASP server predictions](https://drive.google.com/open?id=1BfTn4mImVWqocIs_eh-C4SMUn4siQo72),
- [CASP predictions](https://drive.google.com/open?id=1j3fmwkxWiYHfcYVs1uVzVv9ZEQUVihKw),
- [NMA decoy protein structures](https://drive.google.com/open?id=1a-daGSgKNdHdkr_TkzXaxWE2PD0L7Nzp).

## Get the raw CASP data from the original sources
See [CASP/README](./CASP).
This is a very rich set of protein structure predictions collected over
the years of running the [Critical Assessment of protein Structure Prediction](http://predictioncenter.org/) (CASP) experiments.

### Generate NMADecoys (optional)
One of the ways to enlarge a training set for protein quality assessment is to use the [NOLB tool](https://team.inria.fr/nano-d/software/nolb-normal-modes/).
This is an extremely fast and efficient method for data augmentation, which
proved to increase the accuracy of trained protein quality assessment methods.

See [NMADecoys](./NMADecoys) for more details.

### Compute ground truth for decoy protein models
The [TMscore tool](https://zhanglab.ccmb.med.umich.edu/TM-score/) is used
to measure the similarity between protein structures.
By running it for pairs (protein model, native structure) one can compute ground truth for the quality of the protein models.

*Install:*
```bash
wget https://zhanglab.ccmb.med.umich.edu/TM-score/TMscore.f
gfortran -O3 -ffast-math -lm -o TMscore TMscore.f
```

To run TMscore and obtain formatted scores, run script [compute_scores.py](./compute_scores.py).
```bash
./compute_scores.py <structures_paths> <native_pattern> <decoy_pattern> <num_threads>
```
*Example:*
```bash
./compute_scores.py "CASP/data/CASP10/*" "T0*.pdb" "*" 8
```
