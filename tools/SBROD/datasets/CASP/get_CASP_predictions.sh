#!/usr/bin/bash


wget -r -A SUMMARY.lga_sda_4.txt --no-parent http://predictioncenter.org/download_area/CASP8/results_sda/
mv predictioncenter.org/download_area/CASP8/results_sda/ predictioncenter.org/download_area/CASP8/results_LGA_sda/

wget -r -A SUMMARY.lga_sda_4.txt --no-parent http://predictioncenter.org/download_area/CASP9/results_LGA_sda/

for i in {8..12}; do
    wget -r -A tar.gz --no-parent http://predictioncenter.org/download_area/CASP$i/predictions/
done

rm predictioncenter.org/download_area/CASP*/results_LGA_sda/T*-D*
for x in predictioncenter.org/download_area/CASP*/results_LGA_sda/T*.SUMMARY.lga_sda_4.txt; do
    mv $x ${x%lga_sda_4.txt}lga_sda.txt
done

for i in {10..12}; do
    wget -r -A SUMMARY.lga_sda.txt --no-parent http://predictioncenter.org/download_area/CASP$i/results_LGA_sda/
done

rm -f predictioncenter.org/download_area/CASP10/predictions/TR698.tar.gz

for x in predictioncenter.org/download_area/CASP*/predictions/T*.tar.gz; do
    scores=$(dirname $x)/../results_LGA_sda/$(basename ${x%.tar.gz}).SUMMARY.lga_sda.txt
    if [ -f $scores ]; then
        CASP_dir=data/$(basename $(dirname $(dirname $x)))Predictions
        mkdir -p ${CASP_dir}
        tar -xzf $x -C ${CASP_dir}/ >> errorlog.txt
        cat $scores | sed 's/.lga:SUMMARY(GDT)//g' | sed 's/ \+/\t/g' > ${CASP_dir}/$(basename ${x%.tar.gz})/scores.txt
    fi
done
