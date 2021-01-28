# Smooth Backbone-Reliant Orientation-Dependent Protein Quality Assessment
A single-model protein quality assessment method.

## Reference

> Mikhail Karasikov, Guillaume Pagès, Sergei Grudinin.
Smooth orientation-dependent scoring function for coarse-grained protein quality assessment.
Bioinformatics (2018). [https://doi.org/10.1093/bioinformatics/bty1037](https://doi.org/10.1093/bioinformatics/bty1037).

Developed as a Master project by Mikhail Karasikov (https://github.com/karasikov) in the team of Sergei Grudinin at Inria, Grenoble.

* Does not require modeled side-chains (only backbone)
* Does not depend on the side-chain packing
* Continuous function of coordinates of heavy backbone atoms
* Based only on geometric structural features
* Fast (~1 sec per structure) and competitive to the state-of-the-art

### Download precompiled standalone executables
Precompiled standalone executables are available from here.
* [Linux](https://drive.google.com/open?id=0B3zcrZZIqs3fcE9pSWRWRGpCdmc)
* [MacOS](https://drive.google.com/open?id=0B3zcrZZIqs3fODhGUTU0dFVrTjA)
* [Windows](https://drive.google.com/open?id=0B3zcrZZIqs3fazJtR0JtNjhBd2M)

It is recommended to score multiple structures in one run as intialization of the scoring model in each run takes around 1.8 sec.

Namely, `ls CASP12Stage1/T0859/server* | xargs -n 1 ./sbrod` scores 20 protein models from CASP12 Stage1 (with 133 residues each) separately in 48 sec
while scoring of all protein models in one run with `./sbrod CASP12Stage1/T0859/server*` takes just around 15 sec.

To try out SBROD with a scoring model trained on smaller data, replace scoring model `pipeline.pkl` in `sbrod/base/` with another one downloaded
from [here](https://drive.google.com/open?id=0B3zcrZZIqs3fYVpqZkJpbTUzT2c).
Please use these models only for the purpose of reproducing specific results from Karasikov et al. (2018), and the above links otherwise (in applications and benchmarks).

### Tutorials
See [tutorials](http://w17407.vdi.mipt.ru/protein_scoring/static/manuals.html) to learn the framework of SBROD and to see other examples.


## Training custom scoring models
See [protein_scoring](./protein_scoring) to learn the workflow of SBROD or train your own scoring models.

### Getting the data
See [datasets/README](./datasets/).
Either download [data from original sources](./datasets#download-the-casp-data-from-original-sources) or just download [all data with precomputed scores](./datasets#download-preprocessed-data) at once.

### Training and testing scoring models
Find all performed machine learning experiments in [protein_scoring](./protein_scoring).


## Other
* [Install as a CASP quality assessment server](./server)
* [Compile standalone binary executable](./standalone)

### Assess protein models online
[http://www.karasikov.com/proteins](http://www.karasikov.com/proteins)
