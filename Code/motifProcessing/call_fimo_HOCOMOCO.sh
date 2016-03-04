#!/bin/bash
for i in `seq 1 1 22`
do
    nrOfFimoProcesses=`ps -u david | grep fimo | wc -l`
    while [ 10 -lt $nrOfFimoProcesses ]; do
	echo "sleep cycle"
	sleep 100
	nrOfFimoProcesses=`ps -u david | grep fimo | wc -l`
    done
    nice -n 10 bash ./Code/motifProcessing/run_fimo.sh "Data/motifData/HOCOMOCOv9_AD.meme" $i "HOCOMOCO" &
done
