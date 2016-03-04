#!/bin/bash
fileNames=(`find interimData/motifInstances/HOCOMOCOBEDS | grep ".bed"`)
fileBaseNames=(`find interimData/motifInstances/HOCOMOCOBEDS | grep ".bed" | perl -pe  's/HOCOMOCOBEDS/HOCOMOCOHIGHCONF/' | perl -pe  's/(\.bed)$/_HighConf$1/'`)
i=1
len0=${#fileNames[@]}
len=`expr $len0 - 1`
for i in `seq 0 1 $len`
do
    echo $i
    echo ${fileNames[${i}]}
    more ${fileNames[${i}]} | perl -anF"\t" -e 'if($F[5]<1E-5){print}' > ${fileBaseNames[${i}]}
done