#!/usr/bin/env bash


if [ $# -ne 2 ]; then
  echo "Number of provided input arguments: $#"
  echo "Usage: $0 <num_threads> <dataset>"
  echo "Example: $0 30 \"/dev/shm/CASP/data/CASP*\""
  exit 1
fi

num_processes="$1"

dataset=$(cd "$2"; pwd -P)
mkdir -p binaries
binaries=$(cd binaries; pwd -P)


pushd $(dirname $BASH_SOURCE)

mode="residues"
params="-d 4 -b 10 -a 12 -c 5 -n 0 --skip_errors"
python3.5 featurize_dataset_parallel.py "$dataset/*/*[!t]" $mode "$params" $num_processes
python3.5 -c "from aggregate_dataset import aggregate_dataset; import numpy as np; aggregate_dataset('$dataset', '*$mode$(echo "$params" | tr -d ' ').mat', '$binaries')" &

mode="hbonds"
params="-b 6 -a 6 -c 6 -n 2 --skip_errors"
python3.5 featurize_dataset_parallel.py "$dataset/*/*[!t]" $mode "$params" $num_processes
python3.5 -c "from aggregate_dataset import aggregate_dataset; import numpy as np; aggregate_dataset('$dataset', '*$mode$(echo "$params" | tr -d ' ').mat', '$binaries')" &

mode="solvation"
params="-b 3 -a 2 -c 15 --skip_errors"
python3.5 featurize_dataset_parallel.py "$dataset/*/*[!t]" $mode "$params" $num_processes
python3.5 -c "from aggregate_dataset import aggregate_dataset; import numpy as np; aggregate_dataset('$dataset', '*$mode$(echo "$params" | tr -d ' ').mat', '$binaries')" &

mode="backboneatom"
params="-b 25 -c 7 -n 0 --residue_type_dependent --skip_errors"
python3.5 featurize_dataset_parallel.py "$dataset/*/*[!t]" $mode "$params" $num_processes
python3.5 -c "from aggregate_dataset import aggregate_dataset; import numpy as np; aggregate_dataset('$dataset', '*$mode$(echo "$params" | tr -d ' ').mat', '$binaries')" &

wait

popd
