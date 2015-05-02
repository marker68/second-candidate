/*
 * encoder.h
 *
 *  Created on: 2014/11/20
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef SRC_ENCODER_H_
#define SRC_ENCODER_H_

#include <iostream>
#include <cmath>
#include <cstring>
#include <cassert>
#include "sc_utilities.h"
#include "bucket.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace SC {
/**
 * Encoder class
 * Constrains: kc^mc < 2^31 and N < 2^31, mp<=16, kp <= 256
 * DO NOT BREAK THESE CONSTRAINS!
 */
class Encoder {
protected:
	PQConfig config;
	float * cq, * pq; // code-books
	int size = 0;
	int non_empty_bucket = 0;
	int * L, * pid;
	ushort * cid;
	unsigned char * codes;
	vector<Bucket> ivf;
	int old_mp;
public:
	Encoder();
	virtual ~Encoder();
	void load_codebooks(const char *, const char *, bool);
	void load_encoded_data(const char *, bool);

	template<typename DataType>
	inline void encode(const char *, int, bool);
	template<typename DataType>
	inline void rencode(const char *, int, bool);

	void distribution(bool);

	void output(const char *, const char *, bool);
	void output2(const char *, const char *, bool);
	PQConfig get_config();
	void statistic(bool detail = false);
	void set_mp(int);
	void set_size(size_t);
};

/**
 * An utility to encode a set of vectors for first time use
 * @param filename the location of the raw data file
 * @param offset if each vector has an offset
 * @param verbose to enable verbose mode
 */
template<typename DataType>
inline void Encoder::encode(
		const char * filename,
		int offset,
		bool verbose) {
	// Load data from file
	DataType * data;
	config.N = load_data<DataType>(filename,data,offset,config.dim,verbose);
	if(config.mc <= 0 || config.mp <= 0) {
		cerr << "Nothing to do" << endl;
		return;
	}

	SimpleCluster::init_array(cid, static_cast<size_t>(config.N)
			* static_cast<size_t>(config.mc));
	SimpleCluster::init_array(codes, static_cast<size_t>(config.N)
			* static_cast<size_t>(config.mp));

	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	// Encode data

	double d = 0.0, d_tmp = 0.0;
	int bsc = config.dim/config.mc,
			bsp = config.dim/config.mp;
	int base = 0;
	size_t p = static_cast<size_t>(config.N) / max_threads;

#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for private(base,d,d_tmp)
#endif
		for(int i0 = 0; i0 < max_threads; i0++) {
			DataType * v_tmp;
			float * v_tmp0, * v_tmp1, * v_tmp2, * v_tmp3, * v_tmp4, * v_tmp5;
			if(!SimpleCluster::init_array<float>(v_tmp0,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}
			if(!SimpleCluster::init_array<float>(v_tmp1,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}
			if(!SimpleCluster::init_array<float>(v_tmp2,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}

			// Range definition
			size_t start = p * static_cast<size_t>(i0);
			size_t  end = start + p;
			size_t base_u = start * config.mc;
			size_t base_c = start * config.mp;
			size_t base_p = start;
			size_t i, j, k;
			ushort min_p;

			if(end > config.N || i0 ==  max_threads - 1)
				end = config.N;

			v_tmp = data + start * static_cast<size_t>(config.dim);
			for(i = start; i < end; i++) {
				for(j = 0; j < config.dim; j++) {
					v_tmp0[j] = static_cast<float>(v_tmp[j]);
				}
				v_tmp5 = v_tmp0;

				// Find the coarse quantizer code
				base = 0;
				v_tmp3 = cq; // fixed this pointer to point into the begin of cq
				for(j = 0; j < config.mc; j++) {
					d = DBL_MAX;
					// Linear search
					v_tmp4 = v_tmp3;
					for(k = 0; k < config.kc; k++) {
						d_tmp = SimpleCluster::distance_l2_square<float>(
								v_tmp5, v_tmp4,bsc);
						if(d > d_tmp) {
							d = d_tmp;
							min_p = k;
						}
						v_tmp4 += bsc;
					}
					v_tmp1 = v_tmp3 + min_p * bsc;
					cid[base_u] = min_p;
					// Residual vector
					for(k = 0; k < bsc; k++) {
						v_tmp2[base++] = v_tmp5[k] - v_tmp1[k];
					}
					v_tmp5 += bsc;
					v_tmp3 = v_tmp4;
					base_u++;
				}

				// Encoding
				// Point this to the begin of pq
				v_tmp3 = pq;
				v_tmp4 = v_tmp2;
				for(j = 0; j < config.mp; j++) {
					d = DBL_MAX;
					// Linear search
					for(k = 0; k < config.kp; k++) {
						d_tmp = SimpleCluster::distance_l2_square<float>(
								v_tmp4, v_tmp3,bsp);
						if(d > d_tmp) {
							d = d_tmp;
							min_p = k;
						}
						v_tmp3 += bsp;
					}
					codes[base_c++] = min_p;
					v_tmp4 += bsp;
				}

				v_tmp += config.dim;
			}
		}
#ifdef _OPENMP
	}
#endif
}

/**
 * An utility to re-encode the data vectors.
 * SO you do not need to re-run some steps in the Encoder::encode() method.
 * @param filename the location of the raw data file
 * @param offset if each vector has an offset
 * @param verbose to enable verbose mode
 */
template<typename DataType>
inline void Encoder::rencode(
		const char * filename,
		int offset,
		bool verbose) {
	// Load data from file
	DataType * data;
	config.N = load_data<DataType>(filename,data,offset,config.dim,verbose);
	if(config.mc <= 0 || config.mp <= 0) {
		cerr << "Nothing to do" << endl;
		return;
	}

	SimpleCluster::init_array(codes, static_cast<size_t>(config.N)
			* static_cast<size_t>(config.mp));

	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	// Encode data

	double d = 0.0, d_tmp = 0.0;
	int bsc = config.dim/config.mc,
			bsp = config.dim/config.mp;
	int base = 0;
	size_t p = static_cast<size_t>(config.N) / max_threads;

#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for private(base,d,d_tmp)
#endif
		for(int i0 = 0; i0 < max_threads; i0++) {
			DataType * v_tmp;
			float * v_tmp0, * v_tmp1, * v_tmp2, * v_tmp3, * v_tmp4, * v_tmp5;
			if(!SimpleCluster::init_array<float>(v_tmp0,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}
			if(!SimpleCluster::init_array<float>(v_tmp1,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}
			if(!SimpleCluster::init_array<float>(v_tmp2,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}

			// Range definition
			size_t start = p * static_cast<size_t>(i0);
			size_t  end = start + p;
			size_t base_u = start * config.mc;
			size_t base_c = start * config.mp;
			size_t base_p = start;
			size_t i, j, k;
			ushort min_p;

			if(end > config.N || i0 ==  max_threads - 1)
				end = config.N;

			for(i = start; i < end; i++) {
				v_tmp = data + static_cast<size_t>(config.dim) * pid[base_p++];
				for(j = 0; j < config.dim; j++) {
					v_tmp0[j] = static_cast<float>(v_tmp[j]);
				}
				v_tmp5 = v_tmp0;

				// Find the coarse quantizer code
				base = 0;
				v_tmp3 = cq; // fixed this pointer to point into the begin of cq
				for(j = 0; j < config.mc; j++) {
					min_p = cid[base_u];
					v_tmp1 = v_tmp3 + min_p * bsc;
					// Residual vector
					for(k = 0; k < bsc; k++) {
						v_tmp2[base++] = v_tmp5[k] - v_tmp1[k];
					}
					v_tmp5 += bsc;
					v_tmp3 += config.kc * bsc;
					base_u++;
				}

				// Encoding
				// Point this to the begin of pq
				v_tmp3 = pq;
				v_tmp4 = v_tmp2;
				for(j = 0; j < config.mp; j++) {
					d = DBL_MAX;
					// Linear search
					for(k = 0; k < config.kp; k++) {
						d_tmp = SimpleCluster::distance_l2_square<float>(
								v_tmp4, v_tmp3,bsp);
						if(d > d_tmp) {
							d = d_tmp;
							min_p = k;
						}
						v_tmp3 += bsp;
					}
					codes[base_c++] = min_p;
					v_tmp4 += bsp;
				}
			}
		}
#ifdef _OPENMP
	}
#endif
}
} /* namespace PQLearn */
#endif /* SRC_ENCODER_H_ */
