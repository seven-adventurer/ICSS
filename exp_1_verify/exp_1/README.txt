=============Related Code Run Samples============

*******The code is run in the following strict order:*****************
_0_make_dir/yyy_mkdir.sh: 
           --Submitting shell script to generate two directories:
	./1  This directory stores the results of the experimental runs, mainly containing the results of the propagation range for a certain number of runs   
	./networks This directory stores experimental network data


_1_get_data/yyy_cp_data.sh:
           --Submitting shell script to copy the experimental network to this catalogue (./networks) from ../../_exp_networks/networks_data/*.txt    

_2_setting_up/CLTD_eval:
           --A runnable file which obtains by compiling ... /... /... /c_LTD_new.c

_2_setting_up/yyy_cp_run.sh:
           --Submitting shell script to copy the executable procedure to other directory  (../1/)

_2_setting_up/yyy1_CLTD.sh:
           --Fix some parameters for multiple runs
           --No need to commit this shell script


_3_submit_job/yyy_submit_job.sh:
           --Determine all experimental parameter ranges and submit experiments


*******Description of relevant input parameters:*****************
_2_setting_up/yyy1_CLTD.sh:
	--netid: network number
	--subnetid: Subnetwork number, normally fixed to 1
	--flag_record_st: Flag for whether to record the state of the node being activated (0-No, 1-Yes)
	--flag_rand_p: Flag for whether the edge activation probability is random or not (0-No, 1-Yes)
	--ave_times: Number of simulations of the diffusion process
	--sample_percent: Percentage of nodes sampled
	--cal_times: Number of sampling method runs
	--permu_num_method: Number of permutations of the sampling method, specifically, since there are homo-segments in the node score ordering of the sampling method, we ensure the validity of the method by rearranging the homo-segments
	--active_ratio: Activation probability, i.e. proportion of influential nodes
	--thre_act: Activation thresholds for the LT model
	--method_id: Method number
	--candidate_size: Number of candidate neighbour nodes
	--choose_p: Probability of selecting a node with greater out-degree
	--basicdir: Directory of Experimental Networks

1. method_id (ref: Individual-centralized Seeding Strategy for Influence Maximization in Information-limited Networks):
	-- 1: random
	-- 2: know30HD (to set sample_percent=0.30)
	-- 3: one-hop
	-- 4: one-hopHD
	-- 5: HD (known all nodes information)
	-- 6: the proposed method
	-- 6: the proposed method(beta)