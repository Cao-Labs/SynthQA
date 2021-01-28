# Generate protein decoy structures using the NMA approach

## Data augmentation
Run [generate_decoys.sh](./generate_decoys.sh) to generate NMA decoy protein structures and enlarge your training set. 

*Example:*
```bash
./generate_decoys.sh "../CASP/data/CASP12/T0860/T0*\.pdb" ./CASP12
```

**Prerequisites:**
* NOLB
([standalone executables are available here](https://team.inria.fr/nano-d/software/nolb-normal-modes/))
  - Alexandre Hoffmann & Sergei Grudinin. NOLB : Non-linear rigid block normal mode analysis method.
    Journal of Chemical Theory and Computation, 2017, 13 (5), pp.2123-2134. DOI: 10.1021/acs.jctc.7b00197.

```bash
wget --no-check-certificate https://files.inria.fr/NanoDFiles/Website/Software/NOLB/Linux/NOLB
chmod 755 NOLB
```
Also [available for MacOS](https://files.inria.fr/NanoDFiles/Website/Software/NOLB/MacOS/NOLB).

> ```bash
> ./NOLB native.pdb -r 2.0 -m -s 100
> ```
> This command will generate 100 protein decoys within RMSD up to 2 A with performing post-optimization of the geometry.
>
> 10 basis modes are used by default (controlled with -n N flag)
>
> The generated pdb file with decoy models can be split by the following command
> ```bash
> grep -n 'MODEL\|ENDMDL' native_nlb_decoys.pdb | cut -d: -f 1 | \
> awk '{if(NR%2) printf "sed -n %d,",$1+1; else printf "%dp native_nlb_decoys.pdb > native_nlb_decoy_%03d.pdb\n", $1-1,NR/2;}' |  bash -sf
> ```