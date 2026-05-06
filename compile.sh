#source /opt/intel/bin/compilervars.sh intel64

cd ./source/

make
#make -t

cd ..

mkdir ./bin
mkdir ./output
cp ./source/*mod ./bin/
cp ./source/*out ./bin/
