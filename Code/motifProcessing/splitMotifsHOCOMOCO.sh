#!/bin/bash


for i in `seq 1 1 22`
do
    echo ${i}
    if [ ${i} -eq 1 ];then

	more interimData/motifInstances/all_motifTsHOCOMOCO_${i} | awk -F"\t" '{print > "interimData/motifInstances/HOCOMOCOSPLITTED/"$1"_splitted.txt"}'
    else

	more interimData/motifInstances/all_motifTsHOCOMOCO_${i} | awk -F"\t" '{print >> "interimData/motifInstances/HOCOMOCOSPLITTED/"$1"_splitted.txt"}'
    fi
done
fileNames=(`find interimData/motifInstances/HOCOMOCOSPLITTED/ -name "*.txt"`)
fileNamesOut=(`find interimData/motifInstances/HOCOMOCOSPLITTED/ -name "*.txt" | perl -ple 's/_splitted.txt/.bed/' | perl -ple 's/HOCOMOCOSPLITTED/HOCOMOCOBEDS/'`)
len=`find interimData/motifInstances/HOCOMOCOSPLITTED/ | wc -l`
for i in `seq 1 1 ${len}`
do
    echo ${i}
    more ${fileNames[$i]} | perl -anlF"\t" -e 'print "$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]"' > interimData/tmp
    sort-bed interimData/tmp > ${fileNamesOut[$i]}
done
