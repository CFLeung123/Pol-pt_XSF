#!/bin/bash
######## If using openmp

#Locate the openmp libraries
#source /opt/intel/bin/compilervars.sh intel64

#Number of threads to use:
export OMP_NUM_THREADS=1

#set the stacksize of the job
export OMP_STACKSIZE=2G






#Run the program
for size in `cat sizelist`
#for size in {6..8..2}
 do

echo $size

cat in.dat.src | sed "s/=CFG=/$size/" > ./bin/in.dat

cd ./bin
time ./pt_simple.out > ./output.dat
mv ./output.dat ../output/$size.dat
#time ./pt_simple.out
cd ..

done
