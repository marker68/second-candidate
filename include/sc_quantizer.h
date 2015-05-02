/*
 * sc_quantizer.h
 *
 *  Created on: 2014/12/27
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef SC_QUANTIZER_H_
#define SC_QUANTIZER_H_

#include <iostream>
#include <cstdio>
#include <cassert>
#include <k-means.h>
#include <utilities.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "quantizer.h"
#include "sc_utilities.h"

using namespace std;
using namespace SimpleCluster;

namespace SC
{
template<typename DataType>
class SCQuantizer : public PQQuantizer<DataType>{
public:
	SCQuantizer(int,int,int,int,int,bool);
	virtual ~SCQuantizer();

	/**
	 * Calculate the residual raw data
	 */
	inline void calc_residual_vector(bool);
};

/**
 * Constructors
 * @param filename the path to dataset
 * @param _dim the dimensions number
 * @param _m the number of partitions
 * @param _nsc the number of centers
 * @param verbose enable verbose mode
 */
template<typename DataType>
SCQuantizer<DataType>::SCQuantizer (
		int _dim,
		int _m,
		int _nsc,
		int _os,
		int _k,
		bool verbose)
		: PQQuantizer<DataType>::PQQuantizer(_dim,_m,_nsc,_os,_k,verbose){
}


/**
 * Destructors
 */
template<typename DataType>
SCQuantizer<DataType>::~SCQuantizer () {
}

/**
 * Calculate the residual raw data
 * @param cq the coarse quantizer
 * @param verbose enable verbose mode
 */
template<typename DataType>
inline void SCQuantizer<DataType>::calc_residual_vector(bool verbose) {
	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	int i, j, k0, i0;
	size_t p = static_cast<size_t>(this->N) / static_cast<size_t>(max_threads);

#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for private(i, j, k0,i0)
#endif
		for(i0 = 0; i0 < max_threads; i0++) {
			size_t start = p * static_cast<size_t>(i0);
			size_t end = start + p;
			if(end > this->N || i0 == max_threads - 1)
				end = this->N;
			size_t base = start * static_cast<size_t>(this->dim);
			float * v_tmp1, * v_tmp2, * v_tmp3, * tmp;
			float d, d2, d_tmp;
			int pos, pos2, m, bs, kk, k1;
			for(k0 = 0; k0 < this->k; k0++) {
				tmp = this->data + base;
				m = this->mc[k0];
				bs = this->dim / this->part;
				kk = this->kc[k0];
				for(i = start; i < end; i++) {
					v_tmp1 = this->cq[k0];
					for(j = 0; j < this->part; j++) {
						d = d2 = FLT_MAX;
						v_tmp3 = v_tmp1;
						for(k1 = 0; k1 < kk; k1++) {
							d_tmp = SimpleCluster::distance_l2_square(tmp,v_tmp1,bs);
							if(d > d_tmp) {
								d2 = d;
								d = d_tmp;
								pos2 = pos;
								pos = k1;
							} else {
								if(d2 > d_tmp) {
									d2 = d_tmp;
									pos2 = k1;
								}
							}
							v_tmp1 += bs;
						}

						v_tmp3 += pos2 * bs;
						for(k1 = 0; k1 < bs; k1++) {
							tmp[k1] -= v_tmp3[k1];
						}
						tmp += bs;
					}
				}
			}
		}
#ifdef _OPENMP
	}
#endif
}

} /* namespace PQLearn */

#endif /* SC_QUANTIZER_H_ */
