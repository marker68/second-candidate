/*
 * sc_encoder.h
 *
 *  Created on: 2014/11/20
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef SRC_MREncoder2_H_
#define SRC_MREncoder2_H_

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include "encoder.h"
#include "sc_algorithm.h"
#include "sc_utilities.h"
#include "bucket.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace SC {
/**
 * MREncoder class
 * Constrains: kc^mc < 2^31 and N < 2^31, mp<=16, kp <= 256
 * DO NOT BREAK THESE CONSTRAINS!
 */
class SCEncoder : public Encoder {
protected:
	int nc; // the number of NN to be considered
public:
	SCEncoder() : Encoder::Encoder() {
		nc = 2;
	}
	SCEncoder(int _nc) : Encoder::Encoder() {
		nc = _nc;
		cout << "nc = " << nc << endl;
	}
	virtual ~SCEncoder() {
	}

	template<typename DataType>
	inline void encode(const char *, int, bool);
	template<typename DataType>
	inline void encode2(const char *, int, bool);

	void distribution(bool);

	void output(const char *, const char *, bool);
};

/**
 * An utility to encode a set of vectors for first time use
 * @param filename the location of the raw data file
 * @param offset if each vector has an offset
 * @param verbose to enable verbose mode
 */
template<typename DataType>
inline void SCEncoder::encode2(
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
			* static_cast<size_t>(config.mc << 1));
	SimpleCluster::init_array(codes, static_cast<size_t>(config.N)
			* static_cast<size_t>(config.mp));

	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	// Encode data

	double d = 0.0, d_tmp = 0.0, d2 = 0.0, d_tmp2 = 0.0;
	int bsc = config.dim/config.mc,
			bsp = config.dim/config.mp;
	int base = 0;
	ushort min, min2;
	size_t p = static_cast<size_t>(config.N) / max_threads;

#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for private(base,d,d2,d_tmp,d_tmp2,min,min2)
#endif
		for(int i0 = 0; i0 < max_threads; i0++) {
			DataType * v_tmp;
			float * v_tmp0, * v_tmp1, * v_tmp2, * v_tmp3, * v_tmp4, * v_tmp5;
			ushort * u_tmp1, * u_tmp2;
			unsigned char * c_tmp1;
			if(!SimpleCluster::init_array<float>(v_tmp0,config.dim) ||
					!SimpleCluster::init_array<float>(v_tmp1,config.dim) ||
					!SimpleCluster::init_array<float>(v_tmp2,config.dim)
			) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}

			// Range definition
			size_t start = p * static_cast<size_t>(i0);
			size_t  end = start + p;
			size_t base_u = start * static_cast<size_t>(config.mc << 1);
			u_tmp1 = cid + base_u;
			u_tmp2 = u_tmp1 + config.mc;
			size_t base_c = start * static_cast<size_t>(config.mp);
			c_tmp1 = codes + base_c;
			size_t base_p = start;
			size_t i, j, k;

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
					d = d2 = DBL_MAX;
					// Linear search
					v_tmp4 = v_tmp3;
					for(k = 0; k < config.kc; k++) {
						d_tmp = SimpleCluster::distance_l2_square<float>(
								v_tmp5, v_tmp4,bsc);
						if(d >= d_tmp) {
							d2 = d;
							d = d_tmp;
							min2 = min;
							min = k;
						} else {
							if(d2 >= d_tmp) {
								d2 = d_tmp;
								min2 = k;
							}
						}
						v_tmp4 += bsc;
					}

					*(u_tmp1++) = min;
					*(u_tmp2++) = min2;
					v_tmp1 = v_tmp3 + min * bsc;
					// Residual vector
					for(k = 0; k < bsc; k++) {
						v_tmp2[base++] = v_tmp5[k] - v_tmp1[k];
					}
					v_tmp5 += bsc;
					v_tmp3 = v_tmp4;
				}
				u_tmp1 = u_tmp2;
				u_tmp2 = u_tmp1 + config.mc;

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
							min = k;
						}
						v_tmp3 += bsp;
					}
					*(c_tmp1++) = min;
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
 * An utility to encode a set of vectors for first time use.
 * Use this one if you want to consider more than 2-NNs.
 * @param filename the location of the raw data file
 * @param offset if each vector has an offset
 * @param verbose to enable verbose mode
 */
template<typename DataType>
inline void SCEncoder::encode(
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
			* static_cast<size_t>(config.mc * nc));
	SimpleCluster::init_array(codes, static_cast<size_t>(config.N)
			* static_cast<size_t>(config.mp));

	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	// Encode data
	int bsc = config.dim/config.mc,
			bsp = config.dim/config.mp;
	size_t p = static_cast<size_t>(config.N) / max_threads;

#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for //private()
#endif
		for(int i0 = 0; i0 < max_threads; i0++) {
			float d_tmp[config.kc], d, dt;
			int id[config.kc], min;
			int base = 0;
			DataType * v_tmp;
			float * v_tmp0, * v_tmp1, * v_tmp2, * v_tmp3, * v_tmp4, * v_tmp5;
			ushort ** u_tmp1 = (ushort **)::operator new(config.mc * sizeof(ushort *));
			ushort * u_tmp;
			unsigned char * c_tmp1;
			if(!SimpleCluster::init_array<float>(v_tmp0,config.dim) ||
					!SimpleCluster::init_array<float>(v_tmp1,config.dim) ||
					!SimpleCluster::init_array<float>(v_tmp2,config.dim)
			) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}

			// Range definition
			size_t start = p * static_cast<size_t>(i0);
			size_t  end = start + p;
			size_t base_u = start * static_cast<size_t>(config.mc * nc);
			u_tmp1[0] = cid + base_u;
			size_t i, j, k;
			for(i = 1; i < config.mc; i++)
				u_tmp1[i] = u_tmp1[i-1] + nc;
			size_t base_c = start * static_cast<size_t>(config.mp);
			c_tmp1 = codes + base_c;
			size_t base_p = start;

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
					u_tmp = u_tmp1[j];
					// Linear search
					v_tmp4 = v_tmp3;
					for(k = 0; k < config.kc; k++) {
						d_tmp[k] = SimpleCluster::distance_l2_square<float>(
								v_tmp5, v_tmp4,bsc);
						id[k] = k;
						v_tmp4 += bsc;
					}
					nth_element_id(d_tmp,d_tmp + config.kc,id,nc);
					sort_id(d_tmp,d_tmp+nc,id);
					for(k = 0; k < nc; k++) {
						*(u_tmp++) = id[k];
					}

					v_tmp1 = v_tmp3 + id[0] * bsc;
					// Residual vector
					for(k = 0; k < bsc; k++) {
						v_tmp2[base++] = v_tmp5[k] - v_tmp1[k];
					}
					v_tmp5 += bsc;
					v_tmp3 = v_tmp4;
				}
				u_tmp1[0] += config.mc * nc;
				for(j = 1; j < config.mc; j++)
					u_tmp1[j] = u_tmp1[j-1] + nc;

				// Encoding
				// Point this to the begin of pq
				v_tmp3 = pq;
				v_tmp4 = v_tmp2;
				for(j = 0; j < config.mp; j++) {
					d = FLT_MAX;
					// Linear search
					for(k = 0; k < config.kp; k++) {
						dt = SimpleCluster::distance_l2_square<float>(
								v_tmp4, v_tmp3,bsp);
						if(d > dt) {
							d = dt;
							min = k;
						}
						v_tmp3 += bsp;
					}
					*(c_tmp1++) = min;
					v_tmp4 += bsp;
				}
				v_tmp += config.dim;
			}
		}
#ifdef _OPENMP
	}
#endif
}
} /* namespace PQLearn */
#endif /* SRC_MREncoder_H_ */
