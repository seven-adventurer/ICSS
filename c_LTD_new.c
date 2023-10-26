
#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include "h_rand.h" // large_random_number generator (MT19937, thanks to Makoto Matsumoto and Takuji Nishimura)
#include <time.h>
#include <numeric> 
#include <algorithm>
#include <math.h>
#include <random>
#include <iostream>
using namespace std;



#define FILENAME_LENGHT 1000
#define NMAX 100000000
#define SIRMAX 1000000
// traditional SIR

int nd_idx[NMAX];
char name_a[FILENAME_LENGHT];
int n;
// int re_i_count[SIRMAX],re_s_count[SIRMAX],re_r_count[SIRMAX];
// int all_i_count[SIRMAX],all_s_count[SIRMAX],all_r_count[SIRMAX];
double re_i_count[SIRMAX],re_r_count[SIRMAX];

std::mt19937_64 mtrandm1((unsigned)time(NULL));
std::uniform_real_distribution<double> dicedouble(0.0, 1.00);
std::uniform_int_distribution<int> diceint(1, 1000000000);


int *re_nd_st;

// Output code results as a file
int R_results(int count, int ii, int rep_results, int iden3, int ttt, int infid)
{
    int i;
    FILE *fp_idx;
	
	
    sprintf(name_a,"_%d_%d_%d_%d_%d_Results_i.txt",ii,rep_results, iden3, ttt, infid);
    fp_idx=fopen(name_a,"w+");
    
    fseek(fp_idx, 0, SEEK_SET);
    for(i=0;i<count;i++)
        fprintf(fp_idx, "%f\t", re_i_count[i]);
    fprintf(fp_idx, "\n");
    fclose(fp_idx);
    
    sprintf(name_a,"_%d_%d_%d_%d_%d_Results_r.txt",ii,rep_results, iden3, ttt, infid);
    fp_idx=fopen(name_a,"w+");
    
    fseek(fp_idx, 0, SEEK_SET);
    for(i=0;i<count;i++)
        fprintf(fp_idx, "%f\t", re_r_count[i]);
    fprintf(fp_idx, "\n");
    fclose(fp_idx);
	
	sprintf(name_a,"_%d_%d_%d_%d_%d_Results_nd_st.txt",ii,rep_results, iden3, ttt, infid);
    fp_idx=fopen(name_a,"w+");
    
    fseek(fp_idx, 0, SEEK_SET);
    for(i=0;i<n;i++)
        fprintf(fp_idx, "%d\t", re_nd_st[i]);
    fprintf(fp_idx, "\n");
    fclose(fp_idx);
	
    
    return 0;
}

// Output the results of the code in the form of a file about the propagation range of the selected influence nodes
int R_results_cover(double *re_p, int nn, int ii, int rep_results, int iden3, int flag_rand_p, int iid4, int iid5, int candidate_size, int id_beta, int sample_num)
{
    
    int i;
    FILE *fp_idx;
    sprintf(name_a,"_%d_%d_%d_%d_%d_%d_%d_%d_%d_Results_cover.txt", ii, rep_results, iden3, flag_rand_p, iid4, iid5, candidate_size, id_beta, sample_num);
    fp_idx=fopen(name_a,"w+");
	
    fseek(fp_idx, 0, SEEK_SET);
    for(i=0;i<nn;i++)
        fprintf(fp_idx, "%f\t", re_p[i]);
    fprintf(fp_idx, "\n");
    fclose(fp_idx);
	
	return 0;
}


int R_results_st(int *nd_st, int nn, int ii, int rep_results, int iden3, int flag_rand_p, int iid4, int iid5)
{
    
    int i;
    FILE *fp_idx;
    sprintf(name_a,"__%d_%d_%d_%d_%d_%d_Results_nd_st.txt", ii, rep_results, iden3, flag_rand_p, iid4, iid5);
    fp_idx=fopen(name_a,"w+");
	
    fseek(fp_idx, 0, SEEK_SET);
    for(i=0;i<nn;i++)
        fprintf(fp_idx, "%d\t", nd_st[i]);
    fprintf(fp_idx, "\n");
    fclose(fp_idx);
	
	return 0;
}

int Ini_all(int *nd_st)
{
    return 0;
}

// Count network node out-degrees, and construct inverse (out-degree) networks
int f_network_out(int *adj, int *adj_idx, int *adj_out, int *adj_idx_out, int *DC_out, int nn)
{
	int ii, jj, temp_nei;
	for(ii=0;ii<nn;ii++)
	{
		for(jj=adj_idx[ii]; jj<adj_idx[ii+1]; jj++)
		{
			temp_nei = adj[jj];
			DC_out[temp_nei] += 1;
		}
	}
	adj_idx_out[0] = 0;
	for(ii=0;ii<nn;ii++)
		adj_idx_out[ii+1]=adj_idx_out[ii]+DC_out[ii];

	for(ii=0;ii<nn;ii++)
	{
		for(jj=adj_idx[ii]; jj<adj_idx[ii+1]; jj++)
		{
			
			temp_nei = adj[jj];
			adj_out[adj_idx_out[temp_nei]] = ii;
			adj_idx_out[temp_nei] += 1;
		}
	}
	adj_idx_out[0] = 0;
	for(ii=0;ii<nn;ii++)
		adj_idx_out[ii+1]=adj_idx_out[ii]+DC_out[ii];
	return 1;
}

// random permutation of nodes
void randperm(int *nd_idx, int n)
{
	//srand((unsigned long)time(NULL));
	//sgenrand(rand()+1);
	int i, j, temp_data;
	for(i=0;i<n;i++)
	{
		j = randnew()%(n-i)+i;
		temp_data = nd_idx[j];
		nd_idx[j]=nd_idx[i];
		nd_idx[i]=temp_data;
	}
}

//
void f_get_inmessage_idx(int *adj_old, int *adj_idx_old, int *adj_idx_map, int nn_new)
{
	int ii, jj, zz, temp_nei, flag_temp_check=1;
	for(ii=0;ii<nn_new;ii++)
	{
		for(jj=adj_idx_old[ii]; jj<adj_idx_old[ii+1]; jj++)
		{
			temp_nei=adj_old[jj];
			flag_temp_check=1;
			for(zz=adj_idx_old[temp_nei]; zz<adj_idx_old[temp_nei+1]; zz++)
			{
				if(adj_old[zz]==ii)
				{
					flag_temp_check=0;
					break;
				}
			}
			if(flag_temp_check)
				printf("--------------------error in calculation--------------------------------\n");
			adj_idx_map[jj]=zz;
		}
	}
}

// Initialise the related diffusion arrays
void f_init(int *nd_st, int *act_node, int *inact_node, int *DC_active, int *nd_idx, double *node_act_p, double *node_act_p_2, int nn)
{
	int ii;
	for(ii=0;ii<nn;ii++)
	{
		nd_st[ii]=0;
		act_node[ii]=0;
		inact_node[ii]=nd_idx[ii];
		DC_active[ii]=0;
		node_act_p[ii]=0.0;
		node_act_p_2[ii]=0.0;
	}
}

//Generate normalised random activation probabilities for edges and nodes
void f_generate_p(int *adj, int *adj_idx, double *edge_act_p, double *node_act_p, int nn)
{
	int ii, jj, temp_nei;
	
	for(ii=0;ii<nn;ii++)
		for(jj=adj_idx[ii]; jj<adj_idx[ii+1]; jj++)
		{
			temp_nei=adj[jj];
			edge_act_p[jj]=randnewdouble();
			node_act_p[temp_nei]+=edge_act_p[jj];
		}
	// normalize
	for(ii=0;ii<nn;ii++)
		for(jj=adj_idx[ii]; jj<adj_idx[ii+1]; jj++)
		{
			temp_nei=adj[jj];
			edge_act_p[jj]=edge_act_p[jj]/node_act_p[temp_nei];
		}
}


// Simulate the linear threshold (LT) diffusion model and return the activated nodes
int f_diff(int *adj, int *adj_idx, int *nd_st, int *act_node, int *act_node_new, int *nd_idx, double *edge_act_p, double *node_act_p, double *node_act_p_2, int nn, int act_num, double thre_act)
{
	int ii, jj, old_inact_num=nn-act_num, new_inact_num=nn-act_num, count_level=1, temp_node, temp_nei;
	int act_num_old=act_num, act_num_new, count_act_num=act_num;
	for(ii=0;ii<act_num;ii++)
	{
		temp_node=nd_idx[nn-ii-1];
		nd_st[temp_node]=count_level;
		act_node[ii]=temp_node;
	}
	count_level++;

	while(1)
	{
		act_num_new=0;
		for(ii=0;ii<act_num_old;ii++)
		{
			temp_node=act_node[ii];
			for(jj=adj_idx[temp_node]; jj<adj_idx[temp_node+1]; jj++)
			{
				temp_nei=adj[jj];
				if(nd_st[temp_nei]==0)
				{
					node_act_p[temp_nei]+=edge_act_p[jj];
					// if(node_act_p[temp_nei] > thre_act || fabs(node_act_p[temp_nei] - thre_act)<0.00000001)
					// if(node_act_p[temp_nei] > thre_act)
					if(node_act_p[temp_nei] > thre_act || fabs(node_act_p[temp_nei] - thre_act)<0.00000001)
					{
						// printf("--------node_act_p[temp_nei]: %f-----\n", node_act_p[temp_nei]);
						nd_st[temp_nei]=count_level;
						act_node_new[act_num_new]=temp_nei;
						act_num_new++;
						count_act_num++;
					}
				}
			}
		}

		if(act_num_new==0)
			break;
		for(ii=0;ii<act_num_new; ii++)
			act_node[ii]=act_node_new[ii];
		act_num_old=act_num_new;
		count_level++;
	}
	return count_act_num;
}

// quick sorting algorithm
void quicksort(int *ddate, int *data_idx, int first,int last)
{
	int i, j, pivot, temp;

	if(first<last)
	{
		pivot=first;
		i=first;
		j=last;
		while(i<j)
		{
			while(ddate[i]<=ddate[pivot] && i<last)
				i++;
			while(ddate[j]>ddate[pivot])
				j--;
			if(i<j)
			{
				temp=ddate[i];
				ddate[i]=ddate[j];
				ddate[j]=temp;
				temp=data_idx[i];
				data_idx[i]=data_idx[j];
				data_idx[j]=temp;
			}
		}
		temp=ddate[pivot];
		ddate[pivot]=ddate[j];
		ddate[j]=temp;

		temp=data_idx[pivot];
		data_idx[pivot]=data_idx[j];
		data_idx[j]=temp;
		quicksort(ddate,data_idx,first,j-1);
		quicksort(ddate,data_idx,j+1,last);
	}
}

// permute a partial sequence of nodes
void randperm_part(int *nd_idx, int st_ind, int end_ind)
{
	int i, j, temp_data;
	for(i=st_ind;i<end_ind;i++)
	{
		j = diceint(mtrandm1)%(end_ind-i)+i;
		temp_data = nd_idx[j];
		nd_idx[j]=nd_idx[i];
		nd_idx[i]=temp_data;
	}
}


// 
void f_get_score_frequency(int *nd_score, int *nd_idx, int nn, int act_num)
{
	int ii, jj, thre_score=nd_score[nn-act_num], idx_upper, idx_lower;
	
	for(ii=nn-act_num;ii<nn;ii++)
	{
		if(nd_score[ii]==thre_score)
			idx_upper=ii;
		else
			break;
	}
	
	for(ii=nn-act_num;ii>=0;ii--)
	{
		if(nd_score[ii]==thre_score)
			idx_lower=ii;
		else
			break;
	}
	
	randperm_part(nd_idx, idx_lower, idx_upper+1);
	
	// printf("---------idx_lower: %d, idx_upper: %d, nn-act_num: %d---------\n", idx_lower, idx_upper, nn-act_num);
}


// Map node scores by node number
void f_rank_value(int *temp_nd_score, int *nd_score, int *nd_idx, int nn)
{
	int ii;
	for(ii=0;ii<nn;ii++)
		temp_nd_score[ii]=nd_score[ii];
	
	for(ii=0;ii<nn;ii++)
		nd_score[ii]=temp_nd_score[nd_idx[ii]];
}


// ICSS 
void f_method(int *adj, int *adj_idx, int *nd_idx, int *nd_score, int *DC_in, int *temp_nd_score, int nn, int sample_num, int candidate_size, int *DC_out, int act_num, int *node_state)
{
	int ii, jj, temp_node, temp_nei, temp_idx, temp_dc, max_idx, max_dc, max_nei, temp_st_idx;
	int count_sle_node=0;
	
	for(ii=0;ii<nn;ii++)
	{
		nd_score[ii]=0;
		nd_idx[ii]=ii;
		node_state[ii]=0;
	}

	ii=0;
	while(1)
	{
		ii++;
		temp_node=diceint(mtrandm1)%nn;  // random node
		temp_st_idx=adj_idx[temp_node];

		if(DC_out[temp_node]==0)
			continue;

		max_idx=diceint(mtrandm1)%DC_out[temp_node];  //random neighbor
		max_nei=adj[temp_st_idx+max_idx];
		max_dc=DC_in[max_nei];

		for(jj=1;jj<candidate_size;jj++)  // other neighbor(nums: candidate_size-1)
		{
			temp_idx=diceint(mtrandm1)%DC_out[temp_node];
			temp_nei=adj[temp_st_idx+temp_idx];
			
			temp_dc=DC_in[temp_nei];
			if(temp_dc>max_dc)
			{
				max_dc=temp_dc;
				max_nei=temp_nei;
			}
		}

		nd_score[max_nei]+=1;

		if(node_state[max_nei] == 0)  // non-repeat sampling
		{
			node_state[max_nei]=1;
			count_sle_node+=1;
		}

		if(ii>=sample_num && count_sle_node>=act_num)
			break;
	}
	


	
	// stable_sort(nd_idx, nd_idx+nn, [&nd_score](size_t i1, size_t i2) {return nd_score[i1] < nd_score[i2];});
	// sort the frequency of nodes, and nd_idx
	sort(nd_idx, nd_idx+nn, [&nd_score](size_t i1, size_t i2) {return nd_score[i1] < nd_score[i2];});
	
	// 
	f_rank_value(temp_nd_score, nd_score, nd_idx, nn);
	
	
	// R_results_st(nd_idx, nn, 1, 1, 1, 1, 1, id111);
	
	// R_results_st(nd_score, nn, 101, 101, 101, 101, 101, 101);
	// R_results_st(nd_idx, nn, 111, 111, 111, 111, 111, 111);
	
	// quicksort(nd_score, nd_idx, 0, nn-1);
	
	// R_results_st(nd_idx, nn, 1, 1, 1, 1, 1, id111);
	
	// R_results_st(nd_score, nn, 10, 10, 10, 10, 10, 10);
	// R_results_st(nd_idx, nn, 11, 11, 11, 11, 11, 11);
}

// ICSS_beta 
void f_method_v2(int *adj, int *adj_idx, int *nd_idx, int *nd_score, int *DC_in, int *temp_nd_score, int nn, int sample_num, int candidate_size, int *DC_out, int act_num, int *node_state, double beta)
{
	int ii, jj, temp_node, temp_nei, temp_idx, temp_dc, max_idx, max_dc, max_nei, temp_st_idx;
	int count_sle_node=0;
	
	for(ii=0;ii<nn;ii++)
	{
		nd_score[ii]=0;
		nd_idx[ii]=ii;
		node_state[ii]=0;
	}

	ii=0;
	while(1)
	{
		ii++;
		temp_node=diceint(mtrandm1)%nn;
		temp_st_idx=adj_idx[temp_node];
		if(DC_out[temp_node]==0)
			continue;

		max_idx=diceint(mtrandm1)%DC_out[temp_node];
		max_nei=adj[temp_st_idx+max_idx];
		max_dc=DC_in[max_nei];

		if(dicedouble(mtrandm1)<beta)  // beta probability
		{
			for(jj=1;jj<candidate_size;jj++)
			{
				temp_idx=diceint(mtrandm1)%DC_out[temp_node];
				temp_nei=adj[temp_st_idx+temp_idx];
				
				temp_dc=DC_in[temp_nei];
				if(temp_dc>max_dc)
				{
					max_dc=temp_dc;
					max_nei=temp_nei;
				}
			}
		}

		nd_score[max_nei]+=1;

		if(node_state[max_nei] == 0)
		{
			node_state[max_nei]=1;
			count_sle_node+=1;
		}

		if(ii>=sample_num && count_sle_node>=act_num)
			break;
	}
	
	// sort the frequency of nodes
	sort(nd_idx, nd_idx+nn, [&nd_score](size_t i1, size_t i2) {return nd_score[i1] < nd_score[i2];});
	
	f_rank_value(temp_nd_score, nd_score, nd_idx, nn);
}

// random method
void f_random_strategy(int *nd_idx, int *nd_score, int *temp_nd_score, int nn, int sample_num, int act_num, int *node_state)
{
	int ii, temp_node;
	int count_sle_node=0;
	// randperm_part(nd_idx, 0, nn);
	for(ii=0;ii<nn;ii++)
	{
		nd_score[ii]=0;
		nd_idx[ii]=ii;
		node_state[ii]=0;
	}
	ii=0;
	while(1)
	{
		ii++;
		temp_node=diceint(mtrandm1)%nn;
		nd_score[temp_node]=1;
		if(node_state[temp_node] == 0)
		{
			node_state[temp_node]=1;
			count_sle_node+=1;
		}
		if(ii>=sample_num && count_sle_node>=act_num)
			break;
	}
		
	
	// printf("---------111---------\n");
	// quicksort(nd_score, nd_idx, 0, nn-1);
	// sort(nd_score, nd_score+nn);
	
	sort(nd_idx, nd_idx+nn, [&nd_score](size_t i1, size_t i2) {return nd_score[i1] < nd_score[i2];});
	
	f_rank_value(temp_nd_score, nd_score, nd_idx, nn);
	
}

// know30HD method
void f_know30HD_strategy(int *nd_idx, int *nd_score, int *DC_in, int *temp_nd_score, int nn, int sample_num, int act_num, int *node_state)
{
	int ii, temp_node;
	int count_sle_node=0;
	// randperm_part(nd_idx, 0, nn);
	for(ii=0;ii<nn;ii++)
	{
		nd_score[ii]=0;
		nd_idx[ii]=ii;
		node_state[ii]=0;
	}
	ii=0;
	while(1)
	{
		ii++;
		temp_node=diceint(mtrandm1)%nn;
		nd_score[temp_node]=DC_in[temp_node];
		if(node_state[temp_node] == 0)
		{
			node_state[temp_node]=1;
			count_sle_node+=1;
		}
		if(ii>=sample_num && count_sle_node>=act_num)
			break;
	}

	sort(nd_idx, nd_idx+nn, [&nd_score](size_t i1, size_t i2) {return nd_score[i1] < nd_score[i2];});
	
	f_rank_value(temp_nd_score, nd_score, nd_idx, nn);
}

// one-hop method
void f_one_hop_strategy(int *adj, int *adj_idx, int *nd_idx, int *nd_score, int *DC_out, int *temp_nd_score, int nn, int sample_num, int candidate_size, int act_num, int *node_state)
{
	int ii, jj, temp_node, temp_st_idx, temp_nei, temp_idx;
	int count_sle_node=0;

	for(ii=0;ii<nn;ii++)
	{
		nd_score[ii]=0;
		nd_idx[ii]=ii;
		node_state[ii]=0;
	}

	ii=0;
	while(1)
	{
		ii++;
		temp_node=diceint(mtrandm1)%nn;
		temp_st_idx=adj_idx[temp_node];

		if(DC_out[temp_node]==0)
			continue;

		for(jj=0;jj<candidate_size;jj++)
		{
			temp_idx=diceint(mtrandm1)%DC_out[temp_node];
			temp_nei=adj[temp_st_idx+temp_idx];
			nd_score[temp_nei]+=1;

			if(node_state[temp_nei] == 0)
			{
				node_state[temp_nei]=1;
				count_sle_node+=1;
			}
		}

		if(ii>=sample_num && count_sle_node>=act_num)
			break;
	}
	
	sort(nd_idx, nd_idx+nn, [&nd_score](size_t i1, size_t i2) {return nd_score[i1] < nd_score[i2];});
	
	f_rank_value(temp_nd_score, nd_score, nd_idx, nn);
}

// one-hopHD method
void f_one_hop_HD_strategy(int *adj, int *adj_idx, int *nd_idx, int *nd_score, int *DC_in, int *temp_nd_score, int nn, int sample_num, int candidate_size, int *DC_out, int act_num, int *node_state)
{
	int ii, jj, temp_node, temp_st_idx, temp_nei, temp_idx;
	int count_sle_node=0;
	
	for(ii=0;ii<nn;ii++)
	{
		nd_score[ii]=0;
		nd_idx[ii]=ii;
		node_state[ii]=0;
	}

	ii=0;
	while(1)
	{
		ii++;
		temp_node=diceint(mtrandm1)%nn;
		temp_st_idx=adj_idx[temp_node];

		if(DC_out[temp_node]==0)
			continue;

		for(jj=0;jj<candidate_size;jj++)
		{
			temp_idx=diceint(mtrandm1)%DC_out[temp_node];
			temp_nei=adj[temp_st_idx+temp_idx];
			nd_score[temp_nei]=DC_in[temp_nei];
			if(node_state[temp_nei] == 0)
			{
				node_state[temp_nei]=1;
				count_sle_node+=1;
			}
		}
		if(ii>=sample_num && count_sle_node>=act_num)
			break;
	}

	sort(nd_idx, nd_idx+nn, [&nd_score](size_t i1, size_t i2) {return nd_score[i1] < nd_score[i2];});
	
	f_rank_value(temp_nd_score, nd_score, nd_idx, nn);
}

// HD method
void f_HD_strategy(int *nd_idx, int *nd_score, int *DC_in, int *temp_nd_score, int nn)
{
	int ii, temp_node;
	// randperm_part(nd_idx, 0, nn);
	for(ii=0;ii<nn;ii++)
	{
		nd_score[ii]=DC_in[ii];
		nd_idx[ii]=ii;
	}
	
	sort(nd_idx, nd_idx+nn, [&nd_score](size_t i1, size_t i2) {return nd_score[i1] < nd_score[i2];});
	
	f_rank_value(temp_nd_score, nd_score, nd_idx, nn);
	
}


// Map different methods by method number
void f_method_kitf_method(int *adj, int *adj_idx, int *nd_idx, int *nd_score, int *DC_in, int *DC_out, int *temp_nd_score, int nn, int sample_num,int candidate_size, int method_id, int act_num, int *node_state, double beta)
{
	if(method_id==1)
		f_random_strategy(nd_idx, nd_score, temp_nd_score, nn, sample_num, act_num, node_state);
	else if(method_id==2)
		f_know30HD_strategy(nd_idx, nd_score, DC_in, temp_nd_score, nn, sample_num, act_num, node_state);
	else if(method_id==3)
	{
		// candidate_size=1: one-hop; candidate_size=2: nominate two neighbors per time, and so forth
		// candidate_size=1;
		f_one_hop_strategy(adj, adj_idx, nd_idx, nd_score, DC_out, temp_nd_score, nn, sample_num, candidate_size, act_num, node_state);
	}
	else if(method_id==4)
	{
		// candidate_size=1: one-hop; candidate_size=2: nominate two neighbors per time, and so forth
		candidate_size=1;
		f_one_hop_HD_strategy(adj, adj_idx, nd_idx, nd_score, DC_in, temp_nd_score, nn, sample_num, candidate_size, DC_out, act_num, node_state);
	}
	else if(method_id==5)
		f_HD_strategy(nd_idx, nd_score, DC_in, temp_nd_score, nn);
	else if(method_id==6)
		f_method(adj, adj_idx, nd_idx, nd_score, DC_in, temp_nd_score, nn, sample_num, candidate_size, DC_out, act_num, node_state);
	else if(method_id==7)
		f_method_v2(adj, adj_idx, nd_idx, nd_score, DC_in, temp_nd_score, nn, sample_num, candidate_size, DC_out, act_num, node_state, beta);
	
}


// Repeat many times to estimate the spread of influence nodes
double f_estimate_cover(int *adj, int *adj_idx, int *DC_active, int *nd_st, int *nd_idx, double *node_act_p, double *node_act_p_2, double *edge_act_p, int *act_node, int *inact_node, int nn, int act_num, double thre_act, int flag_rand_p, int ave_times)
{
	int ii, temp_count_act=0;
	double ave_diff_cover=0;
	
	if(flag_rand_p)
	{
		for(ii=0; ii<ave_times; ii++)
		{
			f_init(nd_st, act_node, inact_node, DC_active, nd_idx, node_act_p, node_act_p_2, nn);
			f_generate_p(adj, adj_idx, edge_act_p, node_act_p, nn);
			f_init(nd_st, act_node, inact_node, DC_active, nd_idx, node_act_p, node_act_p_2, nn);
			temp_count_act=f_diff(adj, adj_idx, nd_st, act_node, inact_node, nd_idx, edge_act_p, node_act_p, node_act_p_2, nn, act_num, thre_act);
			ave_diff_cover+=((double)temp_count_act);
		}
		ave_diff_cover=ave_diff_cover/ave_times;
	}
	else
	{
		f_init(nd_st, act_node, inact_node, DC_active, nd_idx, node_act_p, node_act_p_2, nn);
		temp_count_act=f_diff(adj, adj_idx, nd_st, act_node, inact_node, nd_idx, edge_act_p, node_act_p, node_act_p_2, nn, act_num, thre_act);
		ave_diff_cover+=((double)temp_count_act);
	}
	
	return ave_diff_cover;
}



int main(int argc, char *argv[]){
    
    srand((unsigned long)time(NULL));
	
	sgenrand(rand()+1);
	
	
    int ii, jj, t, i, j, k, num_inf, temp_node, sa_idx, ed_idx, temp_nei, t_count, temp_data;
    int *adj;
    int *adj_idx;

    

    int *inf_idx, *new_inf_idx;
	int act_num=1;
    int i_count=0, s_count=0, r_count=0;
    int *nd_st;
    double immu_p, inf_p;
	// Input control parameters
    // double immu_p=atof(argv[4]);
    int immu_n=0;
	int ccount_variable=3;
	
	ccount_variable++;
	int flag_record_st=atoi(argv[ccount_variable]);
	
	ccount_variable++;
	int flag_rand_p=atoi(argv[ccount_variable]);
	
	ccount_variable++;
	int ave_times=atoi(argv[ccount_variable]);
	
	ccount_variable++;
	double sample_percent=atof(argv[ccount_variable]);
	
	ccount_variable++;
	int cal_times=atoi(argv[ccount_variable]);
	
	ccount_variable++;
	int permu_num_method=atoi(argv[ccount_variable]);
	
	ccount_variable++;
	double active_ratio=atof(argv[ccount_variable]);
	
	ccount_variable++;
	double thre_act=atof(argv[ccount_variable]);
	
	ccount_variable++;
    int method_id=atoi(argv[ccount_variable]);

	ccount_variable++;
    double beta=atof(argv[ccount_variable]);
	
	// ccount_variable++;
    // int candidate_size=atoi(argv[ccount_variable]);

	ccount_variable++;
    double candidate_size=atof(argv[ccount_variable]);
	
	ccount_variable++;
    int iden1=atoi(argv[ccount_variable]);
	
	ccount_variable++;
    int iden2=atoi(argv[ccount_variable]);
	
	ccount_variable++;
	int iden3=atoi(argv[ccount_variable]);
    
	
    
    int max_len=0;
    double val_max_len=0.0;
    
    time_t t1,t2;
    time(&t1);
    
	clock_t tStart = clock();
    
	// Read the network structure
    FILE *fp_idx;
    sprintf(name_a, "%s", argv[1]);
	
    fp_idx=fopen(name_a,"rb");

    fseek(fp_idx, 0, SEEK_SET);
    fread(&n,sizeof(int),1,fp_idx);
    fseek(fp_idx, 1*sizeof(int), SEEK_SET);
    adj_idx=(int*)malloc(sizeof(int)*(n+1));

    fread(adj_idx,sizeof(int),n+1,fp_idx);
    fclose(fp_idx);
    
    sprintf(name_a,argv[2]);
    fp_idx=fopen(name_a,"rb");

    fseek(fp_idx, 0, SEEK_SET);
    adj=(int*)malloc(sizeof(int)*adj_idx[n]);
    fread(adj, sizeof(int), adj_idx[n], fp_idx);
    fclose(fp_idx);
    
    sprintf(name_a,argv[3]);
    fp_idx=fopen(name_a,"rb");
    fseek(fp_idx, 0, SEEK_SET);
    fread(nd_idx,sizeof(int),n,fp_idx);
    fclose(fp_idx);
	
	int nn=n, mm=adj_idx[n];
	
    inf_idx=(int*)malloc(sizeof(int)*(n));
	new_inf_idx=(int*)malloc(sizeof(int)*n);
    nd_st=(int*)malloc(sizeof(int)*n);
	re_nd_st=(int*)malloc(sizeof(int)*n);
	
	// ----- new variables -----
	int *DC_in=(int*)malloc(sizeof(int)*n);
	int *DC_active=(int*)malloc(sizeof(int)*n);
	int *adj_idx_map=(int*)malloc(sizeof(int)*mm);
	double *edge_act_p=(double*)malloc(sizeof(double)*mm);
	double *node_act_p=(double*)malloc(sizeof(double)*n);
	double *node_act_p_2=(double*)malloc(sizeof(double)*n);
	int *act_node=(int*)malloc(sizeof(int)*n);
	int *inact_node=(int*)malloc(sizeof(int)*n);
	double *record_cover=(double*)malloc(sizeof(double)*cal_times);
	
	int *nd_score=(int*)malloc(sizeof(int)*n);
	int *temp_nd_score=(int*)malloc(sizeof(int)*n);


	// ----- new variables 2023.03.06 -----
	int *adj_out=(int*)malloc(sizeof(int)*adj_idx[n]);
	int *adj_idx_out=(int*)malloc(sizeof(int)*(n+1));
	int *DC_out=(int*)malloc(sizeof(int)*n);
	int *node_state=(int*)malloc(sizeof(int)*n);
	
	for(ii=0;ii<n;ii++)
	{
		nd_idx[ii]-=1;
		DC_in[ii]=adj_idx[ii+1]-adj_idx[ii];
		DC_out[ii]=0;
		node_state[ii]=0;
	}
	for(ii=0;ii<adj_idx[n];ii++)
		adj[ii]-=1;
	
	f_network_out(adj, adj_idx, adj_out, adj_idx_out, DC_out, nn);




	if(candidate_size<1)
		candidate_size=1;

	
	int flag_sym=1;
	if(flag_rand_p) // Whether the activation probabilities of the edges are random
	{
		int aa=1;
	}
	else
	{

		for(ii=0;ii<nn;ii++)
			for(jj=adj_idx[ii]; jj<adj_idx[ii+1]; jj++)
				edge_act_p[jj] = 1.0 / DC_out[adj[jj]];


	}
	
	if(active_ratio>=1.0)  // Ratio/Number of influential nodes
		act_num=(int)active_ratio;
	else
		act_num=(int)(active_ratio*n);
	
	if(act_num<1)
		act_num=1;


	
	int sample_num;
	if(sample_percent>=1.0)  
		sample_num=(int)(sample_percent);  
	else
		sample_num=(int)(sample_percent*nn);

	int id_beta = (int)(beta*1000);

	
	double ave_diff_cover=0.0, ave_diff_cover_permu=0.0, ave_diff_cover_per=0.0;
	


	int candidate_size_lower = (int)candidate_size, candidate_size_upper=(int)candidate_size+1, curr_candidate_size;  // Integrate the number of candidate sets
	double candidate_p_choose_upper=candidate_size-(double)candidate_size_lower;
	for(ii=0; ii<cal_times; ii++)
	{
		if(dicedouble(mtrandm1)<candidate_p_choose_upper)
			curr_candidate_size=candidate_size_upper;
		else
			curr_candidate_size=candidate_size_lower;

		f_method_kitf_method(adj_out, adj_idx_out, nd_idx, nd_score, DC_in, DC_out, temp_nd_score, nn, sample_num, curr_candidate_size, method_id, act_num, node_state, beta);

		
		ave_diff_cover_permu=0.0;
		for(jj=0;jj<permu_num_method;jj++)
		{
			f_get_score_frequency(nd_score, nd_idx, nn, act_num);
			
			ave_diff_cover_per=f_estimate_cover(adj, adj_idx, DC_active, nd_st, nd_idx, node_act_p, node_act_p_2, edge_act_p, act_node, inact_node, nn, act_num, thre_act, flag_rand_p, ave_times);
			
			ave_diff_cover_permu+=(ave_diff_cover_per/permu_num_method);
		}
		record_cover[ii]=ave_diff_cover_permu;
		ave_diff_cover+=(ave_diff_cover_permu/cal_times);
	}
	
	int iid4=(int)(active_ratio*1000);
	int iid5=(int)(thre_act*1000);
	


	int candidate_size_label=(int)(candidate_size*1000);

	if(active_ratio<1.0)
	{
		if(flag_record_st)
			R_results_st(nd_st, nn, iden1, iden2, iden3, flag_rand_p, iid4, iid5);
		R_results_cover(record_cover, cal_times, iden1, iden2, iden3, flag_rand_p, iid4, iid5, candidate_size_label, id_beta, sample_num);
	}
	else
	{
		if(flag_record_st)
			R_results_st(nd_st, nn, iden1, iden2, iden3, flag_rand_p, act_num, iid5);
		R_results_cover(record_cover, cal_times, iden1, iden2, iden3, flag_rand_p, act_num, iid5, candidate_size_label, id_beta, sample_num);
	}
	
	

	
	
    free(adj_idx);
    adj_idx=NULL;
    free(adj);
    adj=NULL;
    
    free(inf_idx);
    inf_idx=NULL;
    free(new_inf_idx);
    new_inf_idx=NULL;
    free(nd_st);
    nd_st=NULL;
	
	
	free(re_nd_st);
	re_nd_st=NULL;
	
    
	free(DC_in);
	DC_in=NULL;
	free(DC_active);
	DC_active=NULL;
	free(adj_idx_map);
	adj_idx_map=NULL;
	
	free(edge_act_p);
	edge_act_p=NULL;
	free(node_act_p);
	node_act_p=NULL;
	
	free(node_act_p_2);
	node_act_p_2=NULL;
	
	free(act_node);
	act_node=NULL;
	free(inact_node);
	inact_node=NULL;
	
	free(record_cover);
	record_cover=NULL;
	
	free(nd_score);
	nd_score=NULL;
	
	free(temp_nd_score);
	temp_nd_score=NULL;


	free(adj_out);
	adj_out=NULL;
	free(adj_idx_out);
	adj_idx_out=NULL;
	free(DC_out);
	DC_out=NULL;

	free(node_state);
	node_state=NULL;
	
    return 0;
}


