echo "install bigBedToBed"

cd ExternalCode
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod +x ../../ExternalCode/bigBedToBed

echo "install bedops"
if [ `getconf LONG_BIT` -eq "64" ];then
    wget https://github.com/bedops/bedops/releases/download/v2.4.15/bedops_linux_x86_64-v2.4.15.tar.bz2
    tar jxvf bedops_linux_x86_64-v2.4.15.tar.bz2
    rm -r bedops_linux_x86_64-v2.4.15.tar.bz2
fi

if [ `getconf LONG_BIT` -eq "32" ];then
    wget https://github.com/bedops/bedops/releases/download/v2.4.15/bedops_linux_i386-v2.4.15.tar.bz2
    tar jxvf bedops_linux_i386-v2.4.15.tar.bz2
    rm -r bedops_linux_i386-v2.4.15.tar.bz2
fi

wget "http://meme-suite.org/meme-software/4.10.0/meme_4.10.0_2.tar.gz"
##installing meme
tar zxf meme_4.10.0_2.tar.gz
cd meme_4.10.0
myPw=`pwd`
myVar=`echo $myPw | perl -pe 's/ExternalCode\/meme_4\.10\.0//'`
./configure --prefix="${myVar}/ExternalCode/meme" -with-url=http://meme.nbcr.net/meme
make
make test
make install
cd ../
rm meme_4.10.0_2.tar.gz
rm -r meme_4.10.0
