#!/bin/bash

#SBATCH --qos=short
#SBATCH --job-name=xx_job
#SBATCH --account=synet
#SBATCH --array=0-10
#SBATCH --output=xx-%j.out
#SBATCH --error=xx-%j.err
#SBATCH --workdir=./
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

##srun -n 16 hostname

##hostname

net_id=$1

cp_file=yyy1_CLTD.sh
cp_file2=CLTD_eval


target_place=../1/


cp -r ./${cp_file} ${target_place}
cp -r ./${cp_file2} ${target_place}