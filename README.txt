==============Description of relevant documents================
_exp_networks：Store relevant experimental network data
exp_1_verify：For submitting specific experiments and storing the results of experimental runs
c_LTD_new.c: Method Main Code
CLTD_eval：A runnable programme which obtains by compiling c_LTD_new.c
h_rand.c、rand_generator：Functional code for random number generation (not specific, only application)

=========c_LTD_new.c  compile command==========
g++ -std=c++11 -O3 -w c_LTD_new.c rand_generator.c -o CLTD_eval

==========CLTD_eval   run command=============
Reference experiment folder ./exp_1_verify

======CLTD_eval  description of the run results file======
Results file format: _1_1_1_0_10_150_2000_1000_30_Results_cover.txt
_(netid)_(method_id)_(subnetid)_(flag_rand_p)_(int)(active_ratio*1000)_(int)(thre_act*1000)_(candidate_size)_(id_beta)_(sample_num)_Results_cover.txt
