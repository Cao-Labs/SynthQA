#!/bin/bash


function unpack {
  CASP_dataset=$1
  native_structures=$2
  decoy_structures=$3

  mkdir -p data/$CASP_dataset
  for file in predictioncenter.org/download_area/$CASP_dataset/$decoy_structures
  do
    tar -xzf $file -C data/$CASP_dataset/ 2>> errorlog.txt
  done

  if [ $CASP_dataset == 'CASP7' ]; then
    for structure in data/$CASP_dataset/*
    do
      mv $structure data/$CASP_dataset/T$(basename $structure) 2>> errorlog.txt
    done
  fi

  mkdir data/$CASP_dataset/temp
  tar -xzf predictioncenter.org/download_area/$CASP_dataset/$native_structures -C data/$CASP_dataset/temp 2>> errorlog.txt

  if [ $CASP_dataset == 'CASP7' ]; then
    mv data/$CASP_dataset/temp/TARGETS/* data/$CASP_dataset/temp/
    rmdir data/$CASP_dataset/temp/TARGETS
  fi

  rm data/$CASP_dataset/temp/*_* 2>> errorlog.txt
  rm data/$CASP_dataset/temp/T[^0]* 2>> errorlog.txt

  for native in data/$CASP_dataset/temp/*
  do
    native_name=$(basename ${native%.pdb})
    mv $native data/$CASP_dataset/$native_name/$native_name.pdb 2>> errorlog.txt
  done
  rm -r data/$CASP_dataset/temp

  for structure_dir in data/$CASP_dataset/*
  do
    structure=$(basename $structure_dir)
    if [ ! -f $structure_dir/$structure.pdb ]; then
      rm -r $structure_dir
    fi  
  done
}


CASP_dataset='CASP12'
echo $CASP_dataset...
native_structures='targets/casp12.targets_T0.releaseDec022016.tgz'
wget -r -A tar.gz --no-parent http://predictioncenter.org/download_area/$CASP_dataset/server_predictions/
wget -r -A $(basename $native_structures) --no-parent http://predictioncenter.org/download_area/$CASP_dataset/targets/
unpack $CASP_dataset $native_structures "server_predictions/?????.3D.srv.tar.gz"
cp -r predictioncenter.org/download_area/$CASP_dataset/ predictioncenter.org/download_area/${CASP_dataset}Stage1
cp -r predictioncenter.org/download_area/$CASP_dataset/ predictioncenter.org/download_area/${CASP_dataset}Stage2

CASP_dataset='CASP12Stage1'
echo $CASP_dataset...
native_structures='targets/casp12.targets_T0.releaseDec022016.tgz'
unpack $CASP_dataset $native_structures "server_predictions/?????.stage1.3D.srv.tar.gz"

CASP_dataset='CASP12Stage2'
echo $CASP_dataset...
native_structures='targets/casp12.targets_T0.releaseDec022016.tgz'
unpack $CASP_dataset $native_structures "server_predictions/?????.stage2.3D.srv.tar.gz"


CASP_dataset='CASP11'
echo $CASP_dataset...
native_structures='targets/casp11.targets_unsplitted.release11242014.tgz'
wget -r -A tar.gz --no-parent http://predictioncenter.org/download_area/$CASP_dataset/server_predictions/
wget -r -A $(basename $native_structures) --no-parent http://predictioncenter.org/download_area/$CASP_dataset/targets/
unpack $CASP_dataset $native_structures "server_predictions/?????.3D.srv.tar.gz"
cp -r predictioncenter.org/download_area/$CASP_dataset/ predictioncenter.org/download_area/${CASP_dataset}Stage1
cp -r predictioncenter.org/download_area/$CASP_dataset/ predictioncenter.org/download_area/${CASP_dataset}Stage2

CASP_dataset='CASP11Stage1'
echo $CASP_dataset...
native_structures='targets/casp11.targets_unsplitted.release11242014.tgz'
unpack $CASP_dataset $native_structures "server_predictions/?????.stage1.3D.srv.tar.gz"

CASP_dataset='CASP11Stage2'
echo $CASP_dataset...
native_structures='targets/casp11.targets_unsplitted.release11242014.tgz'
unpack $CASP_dataset $native_structures "server_predictions/?????.stage2.3D.srv.tar.gz"


CASP_dataset='CASP10'
echo $CASP_dataset...
native_structures='targets/casp10.targets_unsplitted.noT0695T0739.tgz'
wget -r -A tar.gz --no-parent http://predictioncenter.org/download_area/$CASP_dataset/server_predictions/
wget -r -A $(basename $native_structures) --no-parent http://predictioncenter.org/download_area/$CASP_dataset/targets/
unpack $CASP_dataset $native_structures "server_predictions/?????.3D.srv.tar.gz"
cp -r predictioncenter.org/download_area/$CASP_dataset/ predictioncenter.org/download_area/${CASP_dataset}Stage1
cp -r predictioncenter.org/download_area/$CASP_dataset/ predictioncenter.org/download_area/${CASP_dataset}Stage2

CASP_dataset='CASP10Stage1'
echo $CASP_dataset...
native_structures='targets/casp10.targets_unsplitted.noT0695T0739.tgz'
unpack $CASP_dataset $native_structures "server_predictions/?????.stage1.3D.srv.tar.gz"

CASP_dataset='CASP10Stage2'
echo $CASP_dataset...
native_structures='targets/casp10.targets_unsplitted.noT0695T0739.tgz'
unpack $CASP_dataset $native_structures "server_predictions/?????.stage2.3D.srv.tar.gz"


CASP_dataset='CASP9'
echo $CASP_dataset...
native_structures='targets/casp9.targ_unsplit.tgz'
wget -r -A tar.gz --no-parent http://predictioncenter.org/download_area/$CASP_dataset/server_predictions/
wget -r -A $(basename $native_structures) --no-parent http://predictioncenter.org/download_area/$CASP_dataset/targets/
unpack $CASP_dataset $native_structures "server_predictions/*"


CASP_dataset='CASP8'
echo $CASP_dataset...
native_structures='targets/casp8.targ_unsplit.tar.gz'
wget -r -A tar.gz --no-parent http://predictioncenter.org/download_area/$CASP_dataset/server_predictions/
wget -r -A $(basename $native_structures) --no-parent http://predictioncenter.org/download_area/$CASP_dataset/targets/
unpack $CASP_dataset $native_structures "server_predictions/*"


CASP_dataset='CASP7'
echo $CASP_dataset...
native_structures='targets/targets.all.tgz'
wget -r -A tar.gz --no-parent http://predictioncenter.org/download_area/$CASP_dataset/server_predictions/
wget -r -A $(basename $native_structures) --no-parent http://predictioncenter.org/download_area/$CASP_dataset/targets/
unpack $CASP_dataset $native_structures "server_predictions/*"


CASP_dataset='CASP6'
echo $CASP_dataset...
native_structures='targets_all.tgz'
wget -r -A tar.gz --no-parent http://predictioncenter.org/download_area/$CASP_dataset/MODELS_SUBMITTED/
wget -r -A $native_structures --no-parent http://predictioncenter.org/download_area/$CASP_dataset/
unpack $CASP_dataset $native_structures "MODELS_SUBMITTED/*"
rm data/$CASP_dataset/*/T????[ADRS]*


CASP_dataset='CASP5'
echo $CASP_dataset...
native_structures='casp5_targ.tar.gz'
wget -r -A tar.gz --no-parent http://predictioncenter.org/download_area/$CASP_dataset/MODELS_SUBMITTED/
wget -r -A $native_structures --no-parent http://predictioncenter.org/download_area/$CASP_dataset/
unpack $CASP_dataset $native_structures "MODELS_SUBMITTED/*"
rm data/$CASP_dataset/*/T????[ADRS]*




echo Looking for structures with NaN, it will take some time...
number_of_files=$(find data/*/* | wc -l | tr -d ' ')
i=0

for file in data/*/*/*
do
  if ((i % 17 == 0)); then
    echo -ne "["$(seq -s'=' $((100*i/number_of_files)) | tr -d '[:digit:]')'>'"$(seq -s' ' $((100-100*i/number_of_files)) | tr -d '[:digit:]')""] $i/$number_of_files $(dirname $file)""$(seq -s' ' 10 | tr -d '[:digit:]')"'\r'
  fi
  ((i++))

  if (( $(cat $file | grep nan | wc -l) > 30 )); then
    echo -ne '\r'"$(seq -s' ' 150 | tr -d '[:digit:]')"'\r'
    echo $file:
    cat $file | grep nan | head
    rm $file
    echo $file was removed
  fi
done
echo -ne '\r'"$(seq -s' ' 150 | tr -d '[:digit:]')"'\n'



# remove structures for which TMscore produces errors
rm "data/CASP5/T0137/T0137TS203_1"
rm "data/CASP5/T0149/T0149TS423_1_2"
rm "data/CASP5/T0150/T0150TS326_1"
rm "data/CASP6/T0243/T0243TS094_4"
rm "data/CASP6/T0276/T0276TS532_1"
rm "data/CASP6/T0273/T0273TS406_1"
rm "data/CASP6/T0199/T0199TS331_2"
rm "data/CASP6/T0226/T0226TS506_5"
rm "data/CASP9/T0549/MUSICS_server_TS2"
rm "data/CASP9/T0549/MUSICS_server_TS4"
rm "data/CASP9/T0555/MUSICS_server_TS1"
rm "data/CASP9/T0555/MUSICS_server_TS2"
rm "data/CASP9/T0555/MUSICS_server_TS3"
rm "data/CASP9/T0555/MUSICS_server_TS4"
rm "data/CASP9/T0555/MUSICS_server_TS5"
rm "data/CASP9/T0623/MUSICS-2S_TS1"
rm "data/CASP9/T0623/MUSICS-2S_TS2"
rm "data/CASP9/T0623/MUSICS-2S_TS3"
rm "data/CASP9/T0623/MUSICS-2S_TS4"
rm "data/CASP9/T0623/MUSICS-2S_TS5"
rm "data/CASP9/T0623/MUSICS_server_TS2"
rm "data/CASP9/T0623/MUSICS_server_TS3"
rm "data/CASP9/T0623/MUSICS_server_TS5"
rm "data/CASP9/T0632/MUSICS-2S_TS1"
rm "data/CASP9/T0632/MUSICS-2S_TS2"
rm "data/CASP9/T0632/MUSICS-2S_TS3"
rm "data/CASP9/T0632/MUSICS-2S_TS4"
rm "data/CASP9/T0632/MUSICS-2S_TS5"

