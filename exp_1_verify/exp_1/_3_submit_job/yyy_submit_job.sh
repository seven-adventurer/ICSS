#!/bin/bash

#SBATCH --qos=short
#SBATCH --job-name=xx_job
#SBATCH --account=synet
# #SBATCH --array=0-10
#SBATCH --output=xx-%j.out
#SBATCH --error=xx-%j.err
#SBATCH --workdir=./
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

##srun -n 16 hostname

##hostname


#net_id=$1

# tall_dir=/p/tmp/zhensu/yang_code/c_code/2022.08.05/exps_3/exp_7_CT_comp/
# tall_dir=/p/tmp/zhensu/yang_code/c_code/2022.10.08_IMII/exp_2_verify/exp_1/
# tall_dir=/mnt/e/work_new/vs_code/2022.10.08_IMII/exps/exp_2_verify/exp_1/
tall_dir=../




dddddir=1

sample_percent=0.15

# method_id=1
candidate_size=2

active_ratio_all=(0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10)
thre_act_all=(0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.60)

for net_id in $(seq 5 6)
do
	for method_id in $(seq 8 8)
		do
		for iii in $(seq 0 9)
		do
			active_ratio=${active_ratio_all[${iii}]}
			for jjj in $(seq 0 9)
			do
				thre_act=${thre_act_all[${jjj}]}

				cd ${tall_dir}${dddddir}/
				
				./yyy1_CLTD.sh ${net_id} ${active_ratio} ${thre_act} ${method_id} ${candidate_size} ${sample_percent}
				
				# sbatch yyy1_CLTD.sh ${net_id} ${active_ratio} ${thre_act} ${method_id} ${candidate_size} ${sample_percent}
				
			done
		done
	done
done