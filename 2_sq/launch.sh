#!/bin/bash

for k in 0.1 
do
    
    folderName=k_${k}_$(date +%Y_%m_%d_%H_%M_%S)
    inputFileName=input_k_${k}_$(date +%Y_%m_%d_%H_%M_%S)
    mkdir ${folderName}
    sed "8s/.*/${k}/" input_popsPerturb > ${inputFileName}
    sed -i "10s/.*/${folderName}\//" ${inputFileName}
    
    mpirun -n 4 ./TLGK_LBM_POPS_PERTURB ${inputFileName}
done
