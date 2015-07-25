/*
 * sc_query.h
 *
 *  Created on: 2014/12/19
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef SC_QUERY_H
#define SC_QUERY_H

#include <iostream>
#include <algorithm>
#include <cblas.h>
#include "query.h"

using namespace std;

namespace SC
{

class SCQuery : public PQQuery
{
protected:
	int nc;
public:
	SCQuery();
	SCQuery(int);
	virtual ~SCQuery();
	int load_encoded_data(const char *, bool);

	inline void search_mr_ivf(
			float *,
			float *&, float *&,
			int *, int *&,
			int *&, int *&,int *&, int *&,
			int *&, int *&,
			int *&, int *&, bool *,
			int&, int, int, int, int,
			bool, bool);

	inline void search_mr_ivf3(
			float *,
			float *&, float *&,
			int *, int *&,
			int *&, int *&,int *&, int *&,
			int *&, int *&,
			int *&, //int *&,int *&,
			int *&, int *&, bool *,
			int&, int, int, int, int,
			bool, bool);
};

/**
 * Search method: A demo on single thread mode
 * @param query the query vector
 * @param result the result by identifiers
 * @param R the number of top retrieved results
 * @param verbose to enable verbose mode
 */
inline void SCQuery::search_mr_ivf(float * query,
		float *& v_tmp, float *& dist,
		int * tmp, int *& result,
		int *& hid1, int *& hid2, int *& hid3, int *& hid4,
		int *& s1, int *& s2,
		int *& prebuck, int *& cache, bool * traversed,
		int& sum, int R, int w, int T, int M,
		bool real_dist, bool verbose) {
	if(config.mc != 1 || nc != 2) {
		cerr << "This search method is for MultiRank IVFADC only" << endl;
		return;
	}

	// Temporary pointers: 8 * 8 = 64 bytes
	float * v_tmp1;
	int * i_tmp;
	unsigned char * c_tmp1, * c_tmp2;

	// Step 1: assign the query to coarse quantizer
	size_t i, j, k, l, count = 0, count2 = 0, base = 0, base1, base2, base_c, base_c2, base_k;
	size_t kc = static_cast<size_t>(config.kc);
	size_t kp = static_cast<size_t>(config.kp);
	size_t size2 = kc * kc;
	int h1, h2, h3, h4, m1, m2, c2,
	bid, sum2 = 0, lb = config.kc * config.mc; // 10 * 4 = 40 bytes
	float d_tmp, d_tmp1, d_tmp2, d_tmp3; // 4 bytes
	int bsp = config.dim / config.mp;
	int bsz = bsp * config.kp; // 4 bytes
	int bs = config.kp * config.mp / config.mc;

	pre_compute2(query);
	float q_sum = 0.0;
	v_tmp1 = query;
	for(i = 0; i < config.dim; i++) {
		d_tmp = *(v_tmp1++);
		q_sum += d_tmp * d_tmp;
	}

	v_tmp1 = diff_qc;
	for(j = 0; j < kc; j++) {
		d_tmp = q_sum + (*(v_tmp1++));
		v_tmp[j] = d_tmp;
		tmp[j] = j;
	}

	nth_element_id(v_tmp,v_tmp + config.kc,tmp,M);
	sort_id(v_tmp,v_tmp + M,tmp);

	if(verbose) {
		cout << "Finished STEP 1" << endl;
	}

	// Step 2: Multi-sequences algorithm
	sum = sum2 = count = count2 = 0;
	h3 = tmp[0];
	h4 = tmp[1];
	v_tmp1 = v_tmp + kc;
	d_tmp = v_tmp[0] + v_tmp[1];

	pq_insert4(
			v_tmp1,
			hid1,hid2,hid3,hid4,
			d_tmp,
			0,1,h3,h4,sum2,w,verbose);
	pq_insert4(
			v_tmp1,
			hid1,hid2,hid3,hid4,
			d_tmp,
			1,0,h4,h3,sum2,w,verbose);

	while(count < w && sum < T) {
		if(sum2 > 0) {
			pq_pop4(
					v_tmp1,
					hid1,hid2,hid3,hid4,
					h1,h2,h3,h4,
					sum2,verbose);
			bid = h3 * kc + h4;
			if(bid > 0) {
				l = L[bid] - L[bid-1];
			} else l = 0;
			if(l > 0) {
				prebuck[count] = bid;
				s1[count] = h3;
				//				s2[count] = h4;
				count++;
				sum += l;
			}
			i = kc * h1 + h2;
			traversed[i] = true;
			cache[count2++] = i;
		}


		if(sum < T) {
			j = kc * h2 + h1;
			m1 = i + kc - 1;
			m2 = i - kc + 1;
			if(h1 < kc - 1 &&
					(h2 == 0 || (m1 < size2 && traversed[m1]))
					//					&& h1+1 != h2
			) {
				d_tmp = v_tmp[h1+1] + v_tmp[h2];
				h3 = tmp[h1+1];
				h4 = tmp[h2];

				pq_insert4(
						v_tmp1,
						hid1,hid2,hid3,hid4,
						d_tmp,
						h1+1,h2,h3,h4,sum2,w,verbose);
			}
			if(h2 < kc - 1 &&
					(h1 == 0 || (m2 >= 0 && traversed[m2]))
					//					&& h1 != h2+1
			) {
				d_tmp = v_tmp[h1] + v_tmp[h2+1];
				h3 = tmp[h1];
				h4 = tmp[h2+1];
				pq_insert4(
						v_tmp1,
						hid1,hid2,hid3,hid4,
						d_tmp,
						h1,h2+1,h3,h4,sum2,w,verbose);
			}
		}
	}

	if(verbose) {
		cout << "Finished STEP 2 with " << count << " cells" << endl;
	}

	int count_w = count, count2_w = count2;

	// Step 3: Local search
	// Allocate the memory to store search results
	// Remember to free them after used
	SimpleCluster::init_array(result,sum);
	SimpleCluster::init_array(dist,sum);
	count = 0;

	for(i = 0; i < count_w; i++) {
		bid = prebuck[i];
		h3 = s1[i];
		if(bid > 0) {
			l = L[bid] - L[bid-1];
			i_tmp = pid + L[bid-1];
			c_tmp1 = codes + static_cast<size_t>(L[bid-1]) * config.mp;
		} else {
			l = L[bid];
			i_tmp = pid;
			c_tmp1 = codes;
		}
		if(verbose) {
			cout << "This cell contains " << l << " cadidates with id=" << bid << endl;
		}
		d_tmp = q_sum + diff_qc[h3];
		base1 = h3 * bs;
		// Calculate all l distances
		memcpy(&result[count],i_tmp,l * sizeof(int));
		if(!real_dist) {
			for(j = 0; j < l; j++) {
				base = 0;
				d_tmp1 = d_tmp;
				for(k = 0; k < config.mp; k++) {
					base_c = base + *(c_tmp1++);
					d_tmp1 += (diff_qr[base_c] + dot_cr[base1 + base_c]);
					base += config.kp;
				}
				dist[count++] = d_tmp1;
				//				c_tmp1 += config.mp;
			}
		} else {
			for(j = 0; j < l; j++) {
				dist[count++] = SimpleCluster::distance_l2(
						query,raw_data + i_tmp[j] * config.dim, config.dim);
			}
		}
		if(count >= T) break;
	}

	for(i = 0; i < count2_w; i++) {
		traversed[cache[i]] = 0;
	}

	if(verbose) {
		cout << "Finished STEP 3" << endl;
	}

	// Step 4: Extract the top R
	if(sum >= R) {
		nth_element_id(dist,dist+sum,result,R-1);
		sort_id(dist,dist+R,result);
	}
	if(verbose) {
		cout << "Finished STEP 4" << endl;
	}
}

/**
 * Search method: A demo on single thread mode
 * @param query the query vector
 * @param result the result by identifiers
 * @param R the number of top retrieved results
 * @param verbose to enable verbose mode
 */
inline void SCQuery::search_mr_ivf3(float * query,
		float *& v_tmp, float *& dist,
		int * tmp, int *& result,
		int *& hid1, int *& hid2, int *& hid3, int *& hid4,
		int *& hid5, int *& hid6,
		int *& s1,// int *& s2, int *& s3,
		int *& prebuck, int *& cache, bool * traversed,
		int& sum, int R, int w, int T, int M,
		bool real_dist, bool verbose) {
	if(config.mc != 1 || nc != 3) {
		cerr << "This search method is for MultiRank IVFADC-3 only" << endl;
		return;
	}

	// Temporary pointers: 8 * 8 = 64 bytes
	float * v_tmp1;
	int * i_tmp;
	unsigned char * c_tmp1, * c_tmp2;

	// Step 1: assign the query to coarse quantizer
	size_t i, j, k, l, count = 0, count2 = 0, base = 0, base1, base2, base_c, base_c2, base_k;
	size_t kc = static_cast<size_t>(config.kc);
	size_t kp = static_cast<size_t>(config.kp);
	size_t kc2 = kc * kc;
	size_t size2 = kc * kc * kc;
	int h1, h2, h3, h4, h5, h6, m1, m2, m3,m4,m5,m6,c2,
	bid, sum2 = 0, lb = config.kc * config.mc; // 10 * 4 = 40 bytes
	float d_tmp, d_tmp1, d_tmp2, d_tmp3; // 4 bytes
	int bsp = config.dim / config.mp;
	int bsz = bsp * config.kp; // 4 bytes
	int bs = config.kp * config.mp / config.mc;

	pre_compute2(query);
	float q_sum = 0.0;
	v_tmp1 = query;
	for(i = 0; i < config.dim; i++) {
		d_tmp = *(v_tmp1++);
		q_sum += d_tmp * d_tmp;
	}

	v_tmp1 = diff_qc;
	for(j = 0; j < kc; j++) {
		d_tmp = q_sum + (*(v_tmp1++));
		v_tmp[j] = d_tmp;
		tmp[j] = j;
	}

	nth_element_id(v_tmp,v_tmp + config.kc,tmp,M);
	sort_id(v_tmp,v_tmp + M,tmp);

	if(verbose) {
		cout << "Finished STEP 1" << endl;
	}

	// Step 2: Multi-sequences algorithm
	sum = sum2 = count = count2 = 0;
	h4 = tmp[0];
	v_tmp1 = v_tmp + kc;
	d_tmp = 3.0f * v_tmp[0];// + v_tmp[1] + v_tmp[2];

	pq_insert6(
			v_tmp1,
			hid1,hid2,hid3,hid4,hid5,hid6,
			d_tmp,
			0,0,0,h4,h4,h4,sum2,w,verbose);

	while(count < w && sum < T && count2 < size2) {
		if(sum2 > 0) {
			pq_pop6(
					v_tmp1,
					hid1,hid2,hid3,hid4,hid5,hid6,
					h1,h2,h3,h4,h5,h6,
					sum2,verbose);
			bid = h4 * kc2 + h5 * kc + h6;
			if(bid > 0) {
				l = L[bid] - L[bid-1];
			} else l = 0;
			if(l > 0) {
				prebuck[count] = bid;
				s1[count] = h4;
				count++;
				sum += l;
			}
			i = kc2 * h1 + kc * h2 + h3;
			traversed[i] = true;
			cache[count2++] = i;
		}


		if(sum < T) {
			m1 = i + kc2 - kc;
			m2 = i + kc2 - 1;
			m3 = i + kc - 1;
			m4 = i + kc - kc2;
			m5 = i - kc + 1;
			m6 = i - kc2 + 1;
			if(h1 < kc - 1
					&& (h2 == 0 || (m1 < size2 && traversed[m1]))
					&& (h3 == 0 || (m2 < size2 && traversed[m2]))
			) {
				d_tmp = v_tmp[h1+1] + v_tmp[h2] + v_tmp[h3];
				h4 = tmp[h1+1];
				h5 = tmp[h2];
				h6 = tmp[h3];

				pq_insert6(
						v_tmp1,
						hid1,hid2,hid3,hid4,hid5,hid6,
						d_tmp,
						h1+1,h2,h3,h4,h5,h6,sum2,w,verbose);
			}
			if(h2 < kc - 1
					&& (h3 == 0 || (m3 < size2 && traversed[m3]))
					&& (h1 == 0 || (m4 >= 0 && traversed[m4]))
			) {
				d_tmp = v_tmp[h1] + v_tmp[h2+1] + v_tmp[h3];
				h4 = tmp[h1];
				h5 = tmp[h2+1];
				h6 = tmp[h3];
				pq_insert6(
						v_tmp1,
						hid1,hid2,hid3,hid4,hid5,hid6,
						d_tmp,
						h1,h2+1,h3,h4,h5,h6,sum2,w,verbose);
			}
			if(h3 < kc - 1
					&& (h2 == 0 || (m5 >= 0 && traversed[m5]))
					&& (h1 == 0 || (m6 >= 0 && traversed[m6]))
			) {
				d_tmp = v_tmp[h1] + v_tmp[h2] + v_tmp[h3+1];
				h4 = tmp[h1];
				h5 = tmp[h2];
				h6 = tmp[h3+1];
				pq_insert6(
						v_tmp1,
						hid1,hid2,hid3,hid4,hid5,hid6,
						d_tmp,
						h1,h2,h3+1,h4,h5,h6,sum2,w,verbose);
			}
		}
	}

	if(verbose) {
		cout << "Finished STEP 2 with " << count << " cells" << endl;
	}

	int count_w = count, count2_w = count2;

	// Step 3: Local search
	// Allocate the memory to store search results
	// Remember to free them after used
	SimpleCluster::init_array(result,sum);
	SimpleCluster::init_array(dist,sum);
	count = 0;

	for(i = 0; i < count_w; i++) {
		bid = prebuck[i];
		h3 = s1[i];
		if(bid > 0) {
			l = L[bid] - L[bid-1];
			i_tmp = pid + L[bid-1];
			c_tmp1 = codes + static_cast<size_t>(L[bid-1]) * config.mp;
		} else {
			l = L[bid];
			i_tmp = pid;
			c_tmp1 = codes;
		}
		if(verbose) {
			cout << "This cell contains " << l << " cadidates with id=" << bid << endl;
		}
		d_tmp = q_sum + diff_qc[h3];
		base1 = h3 * bs;
		// Calculate all l distances
		memcpy(&result[count],i_tmp,l * sizeof(int));
		if(!real_dist) {
			for(j = 0; j < l; j++) {
				base = 0;
				d_tmp1 = d_tmp;
				for(k = 0; k < config.mp; k++) {
					base_c = base + *(c_tmp1++);
					d_tmp1 += (diff_qr[base_c] + dot_cr[base1 + base_c]);
					base += config.kp;
				}
				dist[count++] = d_tmp1;
				//				c_tmp1 += config.mp;
			}
		} else {
			for(j = 0; j < l; j++) {
				dist[count++] = SimpleCluster::distance_l2(
						query,raw_data + i_tmp[j] * config.dim, config.dim);
			}
		}
		if(count >= T) break;
	}

	for(i = 0; i < count2_w; i++) {
		traversed[cache[i]] = 0;
	}

	if(verbose) {
		cout << "Finished STEP 3" << endl;
	}

	// Step 4: Extract the top R
	if(sum >= R) {
		nth_element_id(dist,dist+sum,result,R-1);
	}
	if(verbose) {
		cout << "Finished STEP 4" << endl;
	}
}
} /* namespace PQLearn */
#endif /* PQ_MR_QUERY_H_ */
