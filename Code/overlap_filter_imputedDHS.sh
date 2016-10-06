find Data/roadmapDHS
cd Data/roadmapCollected/
#`find ../roadmapDHS | perl -ne '/(E\d+)/ && print "wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-DNase.hotspot.fdr0.01.peaks.bed.gz\n"'`
    `find ../roadmapDHS | perl -ne '/(E\d+)/ && print "wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-DNase.macs2.narrowPeak.gz\n"'`
gunzip *

cat * > allbedresults
sort-bed allbedresults > all_sorted.bed
bedops -m all_sorted.bed > all_merged.bed

cd ../../    

war=(`find Data/roadmapDHS/ |  perl -ne '/(E\d+)/ && print "$1\n"'`)
for el in ${war[@]}
do
echo $el
bedops -e -1 Data/roadmapDHS/${el}-DNase.imputed.narrowPeak.bed Data/roadmapCollected/all_merged.bed > Data/roadmapDHS_filtered/${el}-DNase.imputed.narrowPeak.bed
done

