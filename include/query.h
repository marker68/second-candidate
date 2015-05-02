/*
 * pq_query.h
 *
 *  Created on: 2014/10/15
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef QUERY_H_
#define QUERY_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cstring>
#include <string>
#include <queue>
#include <climits>
#include <cerrno>
#include <cassert>
#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#endif
// Using OpenBLAS for better performance
#include <cblas.h>
#include "sc_utilities.h"
#include "sc_algorithm.h"

using namespace std;

namespace SC {

/*
 * Query class
 * Main jobs are search, update, delete, insert.
 * Search method will be implemented first.
 */
class PQQuery {
protected:
	// Variables
	int size = 0;
	int not_empty = 0;
	float * cq, * pq;
	int * L;
	int * pid;
	unsigned char * codes;
	PQConfig config;
	float * norm_c;
	float * norm_r;
	float * dot_cr;
	float * diff_qr;
	float * diff_qc;
	float * raw_data;
	float * real_dist;
public:
	PQQuery();
	virtual ~PQQuery();
	void load_codebooks(const char *, const char *, bool);
	int load_encoded_data(const char *, bool);
	template<typename DataType>
	inline void load_data(const char *, int, bool);
	inline void pre_compute1();
	inline void pre_compute2(float *);
	inline void pre_compute3(float *);
	inline void search_ivfadc(
			float *,
			float *&, float *&,
			int *&, int *&, int *&,
			int&, int, int,int, bool, bool);
	int get_size();
	int get_full_size();
	double entropy(size_t);
};

template<typename DataType>
inline void PQQuery::load_data(const char * filename, int offset, bool verbose) {
	config.N = SC::load_and_convert_data<DataType,float>(filename,raw_data,offset,config.dim,verbose);
	SimpleCluster::init_array(real_dist,config.N);
}

/**
 * Precompute all things that are query-independent
 */
inline void PQQuery::pre_compute1() {
	SimpleCluster::init_array(norm_c, config.mc * config.kc);
	SimpleCluster::init_array(norm_r, config.mp * config.kp);
	SimpleCluster::init_array(dot_cr, config.kc * config.kp * config.mp);
	SimpleCluster::init_array(diff_qc, config.mc * config.kc);
	SimpleCluster::init_array(diff_qr, config.mp * config.kp);

	float * v_tmp1, * v_tmp2, * v_tmp3;
	float d_tmp, d;
	size_t i, j, i1, j1, k;
	size_t base, base_c, base_p;
	int bsc = config.dim / config.mc;
	int bsp = config.dim / config.mp;

	// Calculate all coarse center norms
	v_tmp1 = cq;
	for(i = 0; i < config.mc * config.kc; i++) {
		d = 0.0;
		for(j = 0; j < bsc; j++) {
			d_tmp = *(v_tmp1++);
			d += d_tmp * d_tmp;
		}
		norm_c[i] = d;
	}

	// Calculate all product center norms
	v_tmp1 = pq;
	for(i = 0; i < config.mp * config.kp; i++) {
		d = 0.0;
		for(j = 0; j < bsp; j++) {
			d_tmp = *(v_tmp1++);
			d += d_tmp * d_tmp;
		}
		norm_r[i] = d;
	}

	// Calculate all dot-products
	base = 0;
	v_tmp1 = cq;
	v_tmp2 = pq;
	for(i = 0; i < config.mc; i++) {
		for(j = 0; j < config.kc; j++) {
			v_tmp3 = v_tmp2;
			for(i1 = 0; i1 < config.mp/config.mc; i1++) {
				for(j1 = 0; j1 < config.kp; j1++) {
					d = 0.0;
					for(k = 0; k < bsp; k++) {
						d += v_tmp1[k] * v_tmp3[k];
					}
					dot_cr[base++] = 2.0f * d;
					v_tmp3 += bsp;
				}
				v_tmp1 += bsp;
			}
		}
		v_tmp2 += config.kp * bsc;
	}
	assert(base == config.kc * config.kp * config.mp);
}

/**
 * Precompute all things that are query dependent
 */
inline void PQQuery::pre_compute2(float * query) {
	size_t i, j, k;
	float *  v_tmp1, * v_tmp2;
	float d, d_tmp = 0.0;
	int bsc = config.dim / config.mc;
	int bsp = config.dim / config.mp;

	v_tmp1 = query;
	v_tmp2 = cq;
	size_t base = 0;
	for(i = 0; i < config.mc; i++) {
		for(j = 0; j < config.kc; j++) {
			d = cblas_sdot(bsc,v_tmp1,1,v_tmp2,1);
			diff_qc[base] = norm_c[base++] - d - d;
			v_tmp2 += bsc;
		}
		v_tmp1 += bsc;
	}

	v_tmp1 = query;
	v_tmp2 = pq;
	base = 0;
	for(i = 0; i < config.mp; i++) {
		for(j = 0; j < config.kp; j++) {
			d = cblas_sdot(bsp,v_tmp1,1,v_tmp2,1);
			diff_qr[base] = norm_r[base++] - d - d;
			v_tmp2 += bsp;
		}
		v_tmp1 += bsp;
	}
}

inline void PQQuery::pre_compute3(float * query) {
	size_t i, j;
	float * v_tmp1 = raw_data;
	float * v_tmp2 = real_dist;
	for(i = 0; i < config.N; i++) {
		*(v_tmp2++) = SimpleCluster::distance_l2_square(query,v_tmp1,config.dim);
		v_tmp1 += config.dim;
	}
}
/**
 * Search method: A demo on single thread mode
 * @param query the query vector
 * @param result the result by identifiers
 * @param R the number of top retrieved results
 * @param verbose to enable verbose mode
 */
inline void PQQuery::search_ivfadc(float * query,
		float *& v_tmp, float *& dist,
		int *& result, int *& buckets, int *& prebuck,
		int& sum, int R, int w, int T, bool real_dist, bool verbose) {
	if(config.mc != 1) {
		cerr << "This search method is for IVFADC only" << endl;
		return;
	}

	// Temporary pointers: 8 * 8 = 64 bytes
	float * v_tmp1;
	float * v_tmp2;
	int * i_tmp;
	unsigned char * c_tmp = codes;

	// Step 1: assign the query to coarse quantizer
	int i, j, k, l, count = 0, count2,
			base = 0, base1, base_c, c, bid, sum2 = 0; // 8 * 4 = 32 bytes
	float d_tmp, d_tmp1; // 4 bytes

	pre_compute2(query);
	float q_sum = 0.0;
	v_tmp1 = query;
	for(i = 0; i < config.dim; i++) {
		d_tmp = *(v_tmp1++);
		q_sum += d_tmp * d_tmp;
	}

	v_tmp1 = diff_qc;
	for(i = 0; i < config.kc; i++) {
		d_tmp = q_sum + *(v_tmp1++);
		pq_insert(v_tmp,buckets,d_tmp,i,sum2,config.kc,verbose);
	}

	if(verbose) {
		cout << "Finished STEP 1" << endl;
		for(i = 0; i < config.kc; i++) {
			cout << buckets[i] << " ";
		}
		cout << endl;
	}


	// Step 2: Local search
	sum2 = config.kc;
	sum = 0;

	for(i = 0; i < w; i++) {
		bid = pq_pop(v_tmp,buckets,sum2,verbose);
		if(bid > 0) {
			l = L[bid] - L[bid-1];
		} else {
			l = L[bid];
		}
		prebuck[count++] = bid;
		sum += l;
		if(sum >= T) break;
	}

	// Allocate the memory to store search results
	// Remember to free them after used
	SimpleCluster::init_array(result,sum);
	SimpleCluster::init_array(dist,sum);
	count = 0;
	int bs = config.dim / config.mp, bsz = bs * config.kp
			,bs1 = config.kp * config.mp / config.mc; // 8 bytes

	for(i = 0; i < w; i++) {
		bid = prebuck[i];
		if(bid > 0) {
			l = L[bid] - L[bid-1];
			i_tmp = pid + L[bid-1];
			c_tmp = codes + static_cast<size_t>(config.mp) *
					static_cast<size_t>(L[bid-1]);
		} else {
			l = L[0];
			i_tmp = pid;
			c_tmp = codes;
		}

		if(verbose) {
			cout << "Searching in bucket " << bid
					<< " that has " << l << " elements" << endl;
		}

		d_tmp1 = q_sum + diff_qc[bid];
		base1 = bid * bs1;

		// Calculate all l distances
		memcpy(&result[count],i_tmp,l * sizeof(int));
//		memcpy(&result[count+sum],i_tmp,l * sizeof(int));
		count2 = count;
		if(!real_dist) {
			for(j = 0; j < l; j++) {
				base = 0;
				d_tmp = d_tmp1;
				for(k = 0; k < config.mp; k++) {
					c = static_cast<int>(*(c_tmp++));
					base_c = base + c;
					d_tmp += (diff_qr[base_c] + dot_cr[base1 + base_c]);
					base += config.kp;
				}
				dist[count] = d_tmp;
				count++;
			}
		} else {
			for(j = 0; j < l; j++) {
				dist[count] = SimpleCluster::distance_l2_square(
						query,raw_data + i_tmp[j] * config.dim,config.dim);
				count++;
			}
		}
//		memcpy(&dist[count2+sum],&dist[count2],l*sizeof(float));
		if(verbose)
			cout << "Searched all " << l << " elements" << endl;
		if(count >= T) break;
	}

	// Step 3: Extract the top R
	if(sum >= R) {
		nth_element_id(dist,dist+sum,result,R-1);
		sort_id(dist,dist+R,result);
	}
}
} /* namespace PQLearn */

#endif /* QUERY_H_ */
