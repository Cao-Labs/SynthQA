#!/usr/bin/env bash

if [ $# -ne 3 ]; then
  echo "Number of provided input arguments: $#"
  echo "Usage: $0 <num_threads> <dataset> <output_path>"
  echo "Example: $0 30 \"/dev/shm/CASP/data/CASP*\""
  exit 1
fi
mkdir -p "$3/binaries"

num_processes="$1"
dataset=$(cd "$2"; pwd -P)
output_path=$(cd "$3"; pwd -P)

binaries=$(cd "$output_path/binaries"; pwd -P)
echo Output path is $output_path


pushd $(dirname $BASH_SOURCE)

mode="residues"
params="-d 4 -b 10 -a 12 -c 5 -n 0"
python3 featurize_dataset_parallel.py "$dataset/*[!t]" $mode "$params" $num_processes "$output_path"
python3 -c "from aggregate_dataset import aggregate_dataset; import numpy as np; aggregate_dataset('$dataset', '*_$mode.mat', '$binaries', '$output_path')" &

mode="hbonds"
params="-b 6 -a 6 -c 6 -n 2"
python3 featurize_dataset_parallel.py "$dataset/*[!t]" $mode "$params" $num_processes "$output_path"
python3 -c "from aggregate_dataset import aggregate_dataset; import numpy as np; aggregate_dataset('$dataset', '*_$mode.mat', '$binaries', '$output_path')" &

mode="solvation"
params="-b 3 -a 2 -c 15"
python3 featurize_dataset_parallel.py "$dataset/*[!t]" $mode "$params" $num_processes "$output_path"
python3 -c "from aggregate_dataset import aggregate_dataset; import numpy as np; aggregate_dataset('$dataset', '*_$mode.mat', '$binaries', '$output_path')" &

mode="backboneatom"
params="-b 25 -c 7 -n 0 --residue_type_dependent"
python3 featurize_dataset_parallel.py "$dataset/*[!t]" $mode "$params" $num_processes "$output_path"
python3 -c "from aggregate_dataset import aggregate_dataset; import numpy as np; aggregate_dataset('$dataset', '*_$mode.mat', '$binaries', '$output_path')" &

wait

popd

