/*
 * multi_query.h
 *
 *  Created on: 2014/10/15
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef MULTI_QUERY_H_
#define MULTI_QUERY_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cstring>
#include <string>
#include <queue>
#include <climits>
#include <cerrno>
#include <cassert>
#include <ctime>
#include <cblas.h>
#include "sc_utilities.h"
#include "sc_algorithm.h"
#include "query.h"

using namespace std;

namespace SC {

/*
 * Query class
 * Main jobs are search, update, delete, insert.
 * Search method will be implemented first.
 */
class MultiQuery : public PQQuery {
protected:
public:
	// The methods of class
	MultiQuery();
	virtual ~MultiQuery();

	inline void search_multi2(
			float *,
			float *&, float *&,float *&,
			int *, int *&,
			int *&, int *&,int *&, int *&,
			int *&, int *&,
			int *&, int *&, bool *,
			int&, int, int,int, int, int&,
			double&, double&,
			bool, bool);
};

/**
 * Search method: A demo on single thread mode
 * @param query the query vector
 * @param result the result by identifiers
 * @param R the number of top retrieved results
 * @param verbose to enable verbose mode
 */
inline void MultiQuery::search_multi2(float * query,
		float *& v_tmp, float *& dist, float *& q,
		int * tmp, int *& result,
		int *& hid1, int *& hid2,int *& hid3, int *& hid4,
		int *& s1, int *& s2,
		int *& prebuck, int *& cache, bool * traversed,
		int& sum, int R, int w, int T, int M, int& e,
		double& t1, double& t2,
		bool real_dist, bool verbose) {
	if(config.mc != 2) {
		cerr << "This search method is for Multi-D-ADC-2 only" << endl;
		return;
	}

	// Temporary pointers: 8 * 8 = 64 bytes
	float * v_tmp1;
	float * v_tmp2;
	float * v_tmp3;
	int * i_tmp;
	int * i_tmp1;
	int * i_tmp2;
	unsigned char * c_tmp;
	clock_t st, ed;

	// Step 1: assign the query to coarse quantizer
	int i, j, k, l, count = 0, base = 0, base1, base2, base_c,
			c = config.mp >> 1,
			h1, h2, h3, h4,
			bid, sum2 = 0, lb = config.kc * config.mc; // 10 * 4 = 40 bytes
	float d_tmp, d_tmp1, d_tmp2; // 4 bytes
	int bsc = config.dim / config.mc;
	int bsp = config.dim / config.mp;
	int bst = bsc * config.kc;
	int bsz = bsp * config.kp; // 4 bytes
	int bs = config.kp * config.mp / config.mc;

	st = clock();
	pre_compute2(query);
	ed = clock();
	t1 += ed - st;
	float q_sum = 0.0;
	v_tmp1 = query;
	for(i = 0; i < config.mc; i++) {
		d_tmp1 = 0.0;
		for(j = 0; j < bsc; j++) {
			d_tmp = *(v_tmp1++);
			d_tmp1 += d_tmp * d_tmp;
		}
		q_sum += d_tmp1;
		q[i] = d_tmp1;
	}

	v_tmp1 = diff_qc;
	for(i = 0; i < config.mc; i++) {
		d_tmp1 = q[i];
		for(j = 0; j < config.kc; j++) {
			d_tmp = d_tmp1 + (*(v_tmp1++));
			v_tmp[base] = d_tmp;
			tmp[base++] = j;
		}
	}

	st = clock();
	nth_element_id(v_tmp,v_tmp + config.kc,tmp,M);
	nth_element_id(v_tmp + config.kc,v_tmp + (config.kc << 1),tmp + config.kc,M);
	sort_id(v_tmp,v_tmp + M,tmp);
	sort_id(v_tmp + config.kc, v_tmp + config.kc + M,tmp + config.kc);
	ed = clock();
	t2 += ed - st;


	// Step 2: Multi-sequences algorithm
	sum = sum2 = count = e = 0;
	int count2 = 0;
	i_tmp1 = tmp;
	i_tmp2 = tmp + config.kc;
	h3 = i_tmp1[0];
	h4 = i_tmp2[0];
	v_tmp1 = v_tmp + lb;
	v_tmp2 = v_tmp;
	v_tmp3 = v_tmp + config.kc;
	d_tmp = v_tmp2[0] + v_tmp3[0];
	pq_insert4(
			v_tmp1,
			hid1,hid2,hid3,hid4,
			d_tmp,0,0,h3,h4,
			sum2,w,verbose);


	while(count2 < w && sum < T) {
		if(sum2 > 0) {
			pq_pop4(
					v_tmp1,
					hid1,hid2,hid3,hid4,
					h1,h2,h3,h4,
					sum2,verbose);
			bid = h3 * config.kc + h4;
			i = h1 * config.kc + h2;
			if(bid > 0) {
				l = L[bid] - L[bid-1];
			} else {
				l = L[bid];
			}
			if(l > 0) {
				prebuck[count] = bid;
				s1[count] = h3;
				s2[count++] = h4;
				sum += l;
			} else e++;
			traversed[i] = true;
			cache[count2++] = i;
		}

		if(sum < T) {
			if(h1 < config.kc - 1 &&
					(h2 == 0 || (i + config.kc <= size && traversed[i + config.kc - 1]))) {
				h3 = i_tmp1[h1+1];
				h4 = i_tmp2[h2];
				d_tmp = v_tmp2[h1+1] + v_tmp3[h2];

				if(i + config.kc < size) {
					pq_insert4(
							v_tmp1,
							hid1,hid2,hid3,hid4,
							d_tmp,h1+1,h2,h3,h4,
							sum2,w,verbose);
				}
			}
			if(h2 < config.kc - 1 &&
					(h1 == 0 || (i + 1 >= config.kc && traversed[i - config.kc + 1]))) {
				h3 = i_tmp1[h1];
				h4 = i_tmp2[h2+1];
				d_tmp = v_tmp2[h1] + v_tmp3[h2+1];

				if(i + 1 < size) {
					pq_insert4(
							v_tmp1,
							hid1,hid2,hid3,hid4,
							d_tmp,h1,h2+1,h3,h4,
							sum2,w,verbose);
				}
			}
		}
	}

	int count_w = count, count2_w = count2;

	// Step 3: Local search
	// Allocate the memory to store search results
	// Remember to free them after used
	SimpleCluster::init_array(result,sum);
	SimpleCluster::init_array(dist,sum);
	count = 0;

	v_tmp1 = diff_qc + config.kc;
	for(i = 0; i < count_w; i++) {
		bid = prebuck[i];
		h3 = s1[i];
		h4 = s2[i];
		if(bid > 0) {
			l = L[bid] - L[bid-1];
			i_tmp = pid + L[bid-1];
			c_tmp = codes + static_cast<size_t>(L[bid-1]) * config.mp;
		} else {
			l = L[bid];
			i_tmp = pid;
			c_tmp = codes;
		}

		d_tmp = q_sum + diff_qc[h3] + v_tmp1[h4];
		base1 = h3 * bs;
		base2 = (h4 + config.kc - 1) * bs;

		// Calculate all l distances
		memcpy(&result[count],i_tmp,l * sizeof(int));
		if(!real_dist) {
			for(j = 0; j < l; j++) {
				base = 0;
				d_tmp1 = d_tmp;
				for(k = 0; k < c; k++) {
					base_c = base + *(c_tmp++);
					d_tmp1 += diff_qr[base_c];
					d_tmp1 += dot_cr[base1 + base_c];
					base += config.kp;
				}
				for(k = c; k < config.mp; k++) {
					base_c = base + *(c_tmp++);
					d_tmp1 += diff_qr[base_c];
					d_tmp1 += dot_cr[base2 + base_c];
					base += config.kp;
				}
				dist[count++] = d_tmp1;
			}
		} else {
			for(j = 0; j < l; j++) {
				dist[count++] = SimpleCluster::distance_l2(
						query,raw_data + i_tmp[j] * config.dim,config.dim);
			}
		}
		if(count >= T) break;
	}

	for(i = 0; i < count2_w; i++) {
		traversed[cache[i]] = 0;
	}

	// Step 3: Extract the top R
	if(sum >= R) {
		nth_element_id(dist,dist+sum,result,R-1);
		sort_id(dist,dist+R,result);
	}
}
} /* namespace PQLearn */

#endif /* MULTI_QUERY_H_ */
