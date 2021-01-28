#!/bin/bash
#sbrod

if [ $# -ne 1 ]; then
  echo "Number of provided input arguments: $#"
  echo "Usage: $0 <num_threads>"
  exit 1
fi

num_processes=$1

dataset="../../datasets/CASP/data/CASP*"
binaries="$(dirname $BASH_SOURCE)/binaries"
mkdir -p $binaries


for dof in "4"; do
  for cutoff in "0" "5" "10" "15" "20"; do
    for d_bins in 6 8 10 12 14; do
      for a_bins in 4 6 9 12 15; do
        for n in 0 1 2; do
          mode="residues"
          params="-d $dof -b $d_bins -a $a_bins -c $cutoff -n $n --skip_errors"

          python3 featurize_dataset_parallel.py "$dataset/*/*[!t]" $mode "$params" $num_processes

          mat_file_pattern="*$mode$(echo "$params" | tr -d ' ').mat"
          python3 -c "from aggregate_dataset import aggregate_dataset; import numpy as np; aggregate_dataset('$dataset', '$mat_file_pattern', '$binaries')" &
        done
      done
    done
  done
done

wait
