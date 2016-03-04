## getting motif data:
cd Data/
#mkdir motifData
cd motifData
wget  "http://autosome.ru/HOCOMOCO/download_helper.php?path=download/HOCOMOCOv9_AD_MEME.txt&name=HOCOMOCOv9_AD_MEME.txt"
mv "download_helper.php?path=download%2FHOCOMOCOv9_AD_MEME.txt&name=HOCOMOCOv9_AD_MEME.txt" HOCOMOCOv9_AD_MEME.txt
more HOCOMOCOv9_AD_MEME.txt | perl -pe 's/\r\n/\n/' > HOCOMOCOv9_AD.meme
cd ../
cd ../
