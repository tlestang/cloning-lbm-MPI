#!/bin/bash

for k in 0.1
do
    #gather nb of proc & Exec name from cla
    NSLOTS=$1
    EXECNAME=$2
    
    
    folderName=k_${k}_$(date +%Y_%m_%d_%H_%M_%S)
    inputFileName=input_k_${k}_$(date +%Y_%m_%d_%H_%M_%S)
    mkdir /scratch/tlestang/${folderName}
    sed "8s/.*/${k}/" input_popsPerturb > ${inputFileName}
    sed -i "10s/.*/${folderName}\//" ${inputFileName}

    cp ${inputFileName} /scratch/tlestang/${folderName}/

    instru_mv_to_scratch="cp ${inputFileName} \/scratch\/tlestang\/"
    instru="mpirun -v -np ${NSLOTS} .\/${EXECNAME} ${inputFileName}"

    batchFileName=batch_k_${k}_$(date +%Y_%m_%d_%H_%M_%S).csh
    jobName="#$ -N ${folderName}"
    nslotsLine="#$ -pe mpi16_debian ${NSLOTS}"
    
    sed "7s/.*/${jobName}/" batch.csh > ${batchFileName}
    sed -i "11s/.*/${nslotsLine}/" ${batchFileName}
    sed -i "38s/.*/${instru_mv_to_scratch}/" ${batchFileName}
    sed -i "45s/.*/${instru}/" ${batchFileName}

    cp ${batchFileName} /scratch/tlestang/${folderName}/

    qsub ${batchFileName}
    
done
