#!/bin/bash

#SBATCH --qos=short
#SBATCH --job-name=xx_job
#SBATCH --account=synet
# #SBATCH --array=3-29
#SBATCH --output=xx-%j.out
#SBATCH --error=xx-%j.err
#SBATCH --workdir=./
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

##srun -n 16 hostname

##hostname

# k=$SLURM_ARRAY_TASK_ID

# k=0

# sleep $k


netid=$1

subnetid=1

flag_record_st=0
flag_rand_p=0
ave_times=20

sample_percent=$6
cal_times=10000
permu_num_method=1

active_ratio=$2
thre_act=$3

method_id=$4
candidate_size=$5
choose_p=$7


./CLTD_eval ${basicdir}net_adj_idx_${netid}.txt ${basicdir}net_adj_${netid}.txt ${basicdir}method_${netid}.txt ${flag_record_st} ${flag_rand_p} ${ave_times} ${sample_percent} ${cal_times}  ${permu_num_method} ${active_ratio} ${thre_act} ${method_id} ${choose_p} ${candidate_size} $netid $method_id $subnetid
basicdir=../networks/


