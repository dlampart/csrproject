mkdir Data/histoneData/
cd Data/histoneData/
mkdir interimData/histoneData/
mkdir interimData/histoneData/rawbeds
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/

ar=(`more index.html | grep href | grep ".narrowPeak.gz" | perl -ne '/(ftp\:.+narrowPeak.gz)">/ && print "$1\n"'`)
for el in ${ar[@]}
do
wget $el
done


for el in `ls | more`
do
echo $el
gunzip $el
done

for el in `ls | more | grep narrowPeak`
do
echo $el
more $el| perl -ne 'if($.>1){print}' > ../../tmp.txt
sort-bed ../../tmp.txt > ../../interimData/histoneData/rawbeds/$el
done
cd  ../../interimData/histoneData/rawbeds
