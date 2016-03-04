#!/bin/bash

memeFile=$1
echo "memeFile:"
echo ${memeFile}
#memeFile="Data/HOCOMOCOv9_AD.meme"
j=$2
infix=$3
echo "chr"
echo ${j}
#j=22

#for j in `seq 22 1 22`
#do
rm interimData/motifInstances/all_motifTs${infix}_$j
DD=(`cat ${memeFile} | perl -ne '/MOTIF/ && print "$.\n"'`)
EE=`expr ${#DD[@]} - 1`
echo ${DD[$EE]}
echo $EE
echo "D[1]: ${DD[1]}"
for i in `seq 0 1 $EE`
do 
echo $i
head -4 ${memeFile} > interimData/motifInstances/myTmpMotifs$j
R=`expr $i`
if [ "$R" -lt "$EE" ]; then
Outer=`expr $i + 1`
echo "Outer: $Outer"
export End=`expr ${DD[$Outer]} - 1`
   else
export End=`more  ${memeFile} | wc -l`
fi
echo "End: $End"
export St=${DD[$i]}
perl -ne 'print if ($.>= "$ENV{'St'}" && $.<= "$ENV{'End'}")' ${memeFile} >> interimData/motifInstances/myTmpMotifs$j
echo "calling"
echo "./ExternalCode/meme/bin/fimo --oc interimData/motifInstances/fimo_out$j --no-qvalue --max-stored-scores 1000000 --thresh 1e-5 interimData/motifInstances/myTmpMotifs$j Data/genomicSequence/chr$j.fa"
./ExternalCode/meme/bin/fimo --oc interimData/motifInstances/fimo_out$j --no-qvalue --max-stored-scores 1000000 --thresh 1e-5 interimData/motifInstances/myTmpMotifs$j Data/genomicSequence/chr$j.fa
cat interimData/motifInstances/fimo_out$j/fimo.txt >> interimData/motifInstances/all_motifTs${infix}_$j
done

