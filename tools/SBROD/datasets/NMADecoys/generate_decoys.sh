#!/bin/bash
# 
# Generate protein decoy structures using the NMA approach
#
# Dependency: NOLB
#   Download standalone executable from
#   https://team.inria.fr/nano-d/software/nolb-normal-modes/
#
# Usage:
#   ./generate_decoys.sh <native_structure_pattern> <output_dir>
#
# Example:
#   ./generate_decoys.sh "../CASP/data/CASP12/*/T0*\.pdb" ./CASP12


native_structures=$1
output_dir=$2


if [ $# -ne 2 ]
then
  echo "Usage: $0 <native_structure_pattern> <output_dir>"
  exit 1
fi


mkdir -p $output_dir


natives=$(ls $native_structures)
echo "Number of native structures:"
echo $(echo $natives | wc -w)


for file in $natives
do
  echo $file
  structure=$(basename $file)
  structure=${structure%.*}
  mkdir -p "$output_dir/$structure"
  cp $file "$output_dir/$structure/$structure.pdb"

  for rmsd in `seq 1.0 9.0` "1.5" "2.5"
  do
    cp "$output_dir/$structure/$structure.pdb" "$output_dir/$structure/${structure}_${rmsd}.pdb"

    ./NOLB "$output_dir/$structure/${structure}_${rmsd}.pdb" -r $rmsd -n 100 -m -s 30 >> log.txt

    decoys_pdb="$output_dir/$structure/${structure}_${rmsd}_nlb_decoys.pdb"
    output_pdb_base="$output_dir/$structure/${structure}_${rmsd}_nlb_decoy"

    grep -n 'MODEL\|ENDMDL' $decoys_pdb | cut -d: -f 1 | \
          awk '{if(NR%2) printf "sed -n %d,",$1+1; else printf "%dp '"'${decoys_pdb}'"' > '"'${output_pdb_base}'"'_%03d.pdb\n", $1-1,NR/2;}' |  bash -sf
    rm $decoys_pdb

    rm "$output_dir/$structure/${structure}_${rmsd}.pdb"
  done
done
