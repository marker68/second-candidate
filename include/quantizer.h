/*
 * quantizer.h
 *
 *  Created on: 2014/06/27
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef QUANTIZER_H_
#define QUANTIZER_H_

#include <iostream>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <unordered_map>
#include <cstring>
#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <k-means.h>
#include <kd-tree.h>
#include "sc_utilities.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace SimpleCluster;

/**
 * Main namespace for Product Quantization implementation
 */
namespace SC {

/**
 * A PQ quantizer class
 *
 * m=1 will create a coarse quantizer
 */
template<typename DataType>
class PQQuantizer {
protected:
	int dim, part, nsc, offset;
	int N; // the number of points data
	float * centers; // size: nsc * dim
	int * labels; // for uses of k-means++
	unsigned int * ulabels;
	float ** cq;
	int k;
	int * kc, * mc;
	int type = 1; // 1: VLFeat; 2: Simple-Cluster

public:
	float * data; // raw vector data; size: N * dim
	double error = 0.0;

	/**
	 * Constructors and destructors
	 */
	PQQuantizer(int,int,int,int,bool);
	PQQuantizer(int,int,int,int,int,bool);
	virtual ~PQQuantizer ();

	/**
	 * Load raw data from another dataset
	 */
	inline void load_data(const char *, bool);
	inline void load_data_mat(float *, int);

	/**
	 * Create sub-quantizers
	 * Using k-means++
	 */
	inline void create_sub_quantizers(bool);

	/**
	 * Distortion
	 */
	inline double distortion(bool);

	/**
	 * Size of raw data
	 */
	inline int size();

	/**
	 * Getters and setters
	 */
	inline float * get_center_at(int, bool);
	inline int * get_label_at(int, bool);
	inline int get_type();
	inline void set_params(int,int,int);
	inline void set_type(int);

	/**
	 * Common output method
	 */
	inline void output(const char *, bool);

	inline void load_codebooks(char **, bool);

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
PQQuantizer<DataType>::PQQuantizer (
		int _dim,
		int _m,
		int _nsc,
		int _os,
		bool verbose) {
	dim = _dim;
	data = nullptr;
	offset = _os;
	N = k = 0;
	part = _m;
	nsc = _nsc;
	labels = nullptr;
	ulabels = nullptr;
	cq = nullptr;
	kc = nullptr;
	mc = nullptr;
	if(!SimpleCluster::init_array<float>(centers,nsc*dim)) {
		if(verbose)
			cerr << "There are some errors occurred while initializing data" << endl;
		exit(EXIT_FAILURE);
	}
	if(dim <= 0 || part <= 0 || dim % part != 0) {
		if(verbose)
			cerr << "Error at parameters of PQQuantizer: dim = " << dim << "; part = " << part << endl;
		exit(EXIT_FAILURE);
	}
}

/**
 * Constructors
 * @param filename the path to dataset
 * @param _dim the dimensions number
 * @param _m the number of partitions
 * @param _nsc the number of centers
 * @param verbose enable verbose mode
 */
template<typename DataType>
PQQuantizer<DataType>::PQQuantizer (
		int _dim,
		int _m,
		int _nsc,
		int _os,
		int _k,
		bool verbose)
		: PQQuantizer<DataType>::PQQuantizer(_dim,_m,_nsc,_os,verbose){
	k = _k;
	cq = (float **)::operator new(k * sizeof(float *));
	init_array(kc,k);
	init_array(mc,k);
}
/**
 * Destructors
 */
template<typename DataType>
PQQuantizer<DataType>::~PQQuantizer () {
	::delete centers;
	centers = nullptr;
	::delete labels;
	labels = nullptr;
	::delete ulabels;
	ulabels = nullptr;
	::delete data;
	data = nullptr;
	for(int i = 0; i < k; i++) {
		::delete cq[i];
		cq[i] = nullptr;
	}
	::delete cq;
	cq = nullptr;
	::delete kc;
	kc = nullptr;
	::delete mc;
	mc = nullptr;
}

/**
 * Load data from disk
 */
template<typename DataType>
inline void PQQuantizer<DataType>::load_data(
		const char * filename,
		bool verbose) {
	N = load_and_convert_data<DataType,float>(filename,data,offset,dim,verbose);
	if(type == 1) {
		if(!SimpleCluster::init_array<unsigned int>(ulabels,N*part)) {
			if(verbose)
				cerr << "There are some errors occurred while initializing data" << endl;
			exit(EXIT_FAILURE);
		}
	} else {
		if(!SimpleCluster::init_array<int>(labels,N*part)) {
			if(verbose)
				cerr << "There are some errors occurred while initializing data" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if (verbose)
		cout << "Initialized " << N << " vectors" << endl;
}

/**
 * Load raw data from another dataset
 * @param _data the data to be copied
 * @param _N the number of vectors
 * @return nothing
 */
template<typename DataType>
inline void PQQuantizer<DataType>::load_data_mat(
		float * _data,
		int _N) {
	if(_N > 0) N = _N;
	data = _data;
}

/**
 * Get size of raw data
 */
template<typename DataType>
inline int PQQuantizer<DataType>::size() {
	return N;
}

/**
 * Sub-quantizers
 * This function will create codes, centers and labels
 * @param verbose enable verbose mode
 */
template<typename DataType>
inline void PQQuantizer<DataType>::create_sub_quantizers(
		bool verbose) {
	int bs = dim / part; //block size
	float * _centers, * _seeds = nullptr;
	int * _labels;
	unsigned int * _ulabels;
	float * tmp_data, *  distances;
	double energy;
	int i, j, col_st;
	KmeansCriteria criteria = {2.0,1.0,1000};
	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	_centers = centers;

	if(type == 1)
		_ulabels = ulabels;
	else
		_labels = labels;

	int start = 0, end = dim / part;
	float d = 0.0f, tmp;
	size_t t;
	for(i = 0; i < part; i++) {
		strip_matrix<float>(data,tmp_data,N,dim,start,end,verbose);
		if(verbose)
			cout << "Creating codebook " << i  << "/" << part << endl;
		if(type == 1) {
			vl_kmeans_exec(
					tmp_data, _centers, _ulabels, distances,
					N, nsc, bs, 1, VlDistanceL2, energy, verbose);
			// Compute the distortion
			for(t = 0; t < N; t++) {
				tmp = distances[t];
				d += tmp;
			}
		} else {
			greg_kmeans<float>(
					tmp_data,_centers,_labels,_seeds,
					KmeansType::KMEANS_PLUS_SEEDS,
					criteria,
					DistanceType::NORM_L2,
					EmptyActs::SINGLETON,
					N,nsc,bs,max_threads,
					verbose);
		}
		start = end;
		end += dim / part;
		if(end > dim) end = dim;
		_centers += nsc * dim / part;
		if(type == 1)
			_ulabels += N;
		else
			_labels += N;
		if(verbose)
			cout << "Finished subcodebook " << i << endl;
		::delete tmp_data;
		tmp_data = nullptr;
	}

	error = sqrt(d);
	return;
}

/**
 * Calculate the distortion of the quantization
 * @param verbose enable verbose mode
 */
template<typename DataType>
inline double PQQuantizer<DataType>::distortion(bool verbose) {
	if(type == 1) return DBL_MAX;
	double e = 0.0, e_tmp;
	float * tmp_data;
	int i;
	float * _centers = centers;
	int * _labels = labels;
	int start = 0, end = dim / part;
	for(i = 0; i < part; i++) {
		strip_matrix(data,tmp_data,N,dim,start,end,verbose);
		e_tmp = SimpleCluster::distortion<float>(tmp_data,_centers,_labels,
				DistanceType::NORM_L2,dim/part,N,nsc,verbose);
		e += e_tmp * e_tmp;
		start = end;
		end += dim / part;
		if(end > dim) end = dim;
		_centers += nsc * dim / part;
		_labels += N;
		::delete tmp_data;
		tmp_data = nullptr;
	}
	error = sqrt(e);
	return sqrt(e);
}

/**
 * Get the centers of a partition data
 * @param id index of the partition
 * @param verbose enable verbose mode
 * @return the centers set of partition id
 */
template<typename DataType>
inline float * PQQuantizer<DataType>::get_center_at(int id, bool verbose) {
	return centers + static_cast<size_t>(id)
			* static_cast<size_t>(nsc) * static_cast<size_t>(dim / part);
}

/**
 * Get all labels of a partition data
 * @param id index of the partition
 * @param verbose enable verbose mode
 * @return the labels set of partition id
 */
template<typename DataType>
inline int * PQQuantizer<DataType>::get_label_at(int id, bool verbose) {
	if(type == 1)
		return ulabels +  id * N;
	return labels + id * N;
}

template<typename DataType>
inline void PQQuantizer<DataType>::set_params(
		int _k,
		int _m,
		int _nsc) {
	for(int i = 0; i < k; i++) {
		::delete cq[i];
		cq[i] = nullptr;
	}
	::delete cq;
	cq = nullptr;
	::delete kc;
	kc = nullptr;
	::delete mc;
	mc = nullptr;
	k = _k;
	cq = (float **)::operator new(k * sizeof(float *));
	SimpleCluster::init_array(kc,k);
	SimpleCluster::init_array(mc,k);
	part = _m;
	nsc = _nsc;
	::delete centers;
	centers = nullptr;
	::delete labels;
	labels = nullptr;
	::delete ulabels;
	ulabels = nullptr;
	if(!SimpleCluster::init_array<float>(centers,nsc*dim)) {
		cerr << "There are some errors occurred while initializing data" << endl;
		exit(EXIT_FAILURE);
	}
	if(dim <= 0 || part <= 0 || dim % part != 0) {
		cerr << "Error at parameters of PQQuantizer" << endl;
		exit(EXIT_FAILURE);
	}
	if(type == 1) {
		if(!SimpleCluster::init_array<unsigned int>(ulabels,N*part)) {
			cerr << "There are some errors occurred while initializing data" << endl;
			exit(EXIT_FAILURE);
		}
	} else {
		if(!SimpleCluster::init_array<int>(labels,N*part)) {
			cerr << "There are some errors occurred while initializing data" << endl;
			exit(EXIT_FAILURE);
		}
	}
}

template<typename DataType>
inline void PQQuantizer<DataType>::set_type(
		int _type) {
	type = _type;
}

template<typename DataType>
inline int PQQuantizer<DataType>::get_type() {
	return type;
}

/**
 * Output the codebook's centers to binary file
 * @param filename the location of output file
 * @param verbose enable verbose mode
 */
template<typename DataType>
inline void PQQuantizer<DataType>::output(
		const char * filename,
		bool verbose) {
	char fn[256];
	// First, write the centers to file named filename.ctr_
	sprintf(fn,"%s.ctr_",filename);
#ifdef _WIN32
	ofstream output;
	output.open(fn, ios::out | ios::binary);
	output << nsc << part << dim;
	int i;
	int base = 3;
	float * tmp = centers;
	for(int i = 0; i < dim * nsc; i++) {
		output << *tmp;
		tmp++;
	}
	output.close();
#else
	int fd = open(fn, O_RDWR | O_CREAT | O_TRUNC,(mode_t)0600); // file description
	if(fd < 0) {
		cerr << "Cannot open the file " << fn << endl;
		exit(EXIT_FAILURE);
	}

	size_t size = (dim * nsc + 3) * sizeof(float);
	int result = lseek(fd, size, SEEK_SET);
	if (result == -1) {
		close(fd);
		cerr << "Error calling lseek() to 'stretch' the file" <<  endl;
		exit(EXIT_FAILURE);
	}

	int status = write(fd, "", 1);
	if(status != 1) {
		cerr << "Cannot write to file" << endl;
		exit(EXIT_FAILURE);
	}

	float * fd_map = (float *)mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (fd_map == MAP_FAILED) {
		close(fd);
		cerr << "Error mmapping the file" << endl;
		exit(EXIT_FAILURE);
	}
	fd_map[0] = nsc;
	fd_map[1] = part;
	fd_map[2] = dim;
	int base = 3;
	float * tmp = centers;

	for(int i = 0; i < dim * nsc; i++) {
		fd_map[base++] = *tmp;
		tmp++;
	}

	if (munmap(fd_map, size) == -1) {
		cerr << "Error un-mmapping the file" << endl;
	}
	close(fd);
#endif
}

/**
 * Load codebooks from disk
 */
template<typename DataType>
inline void PQQuantizer<DataType>::load_codebooks(
		char ** filename,
		bool verbose) {
	// Load all precomputed codebooks into memory
	int i;
	PQConfig config;
	for(i = 0; i < k; i++) {
		load_codebook<float>(filename[i],config,cq[i],0,verbose);
		kc[i] = config.kc;
		mc[i] = config.mc;
	}
}

/**
 * Calculate the residual raw data
 * @param cq the coarse quantizer
 * @param verbose enable verbose mode
 */
template<typename DataType>
inline void PQQuantizer<DataType>::calc_residual_vector(bool verbose) {
	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	int i, j, k0, i0;
	size_t p = static_cast<size_t>(N) / static_cast<size_t>(max_threads);

#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for private(i, j, k0,i0)
#endif
		for(i0 = 0; i0 < max_threads; i0++) {
			size_t start = p * static_cast<size_t>(i0);
			size_t end = start + p;
			if(end > N || i0 == max_threads - 1)
				end = N;
			size_t base = start * static_cast<size_t>(dim);
			float * v_tmp1, * v_tmp2, * tmp;
			float d, d_tmp;
			int pos, m, bs, kk, k1;
			for(k0 = 0; k0 < k; k0++) {
				tmp = data + base;
				m = mc[k0];
				bs = dim / m;
				kk = kc[k0];
				for(i = start; i < end; i++) {
					v_tmp1 = cq[k0];
					for(j = 0; j < m; j++) {
						d = FLT_MAX;
						v_tmp2 = v_tmp1;
						for(k1 = 0; k1 < kk; k1++) {
							d_tmp = SimpleCluster::distance_l2(tmp,v_tmp1,bs);
							if(d > d_tmp) {
								d = d_tmp;
								pos = k1;
							}
							v_tmp1 += bs;
						}

						v_tmp2 += pos * bs;
						for(k1 = 0; k1 < bs; k1++) {
							tmp[k1] -= v_tmp2[k1];
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
}

#endif /* QUANTIZER_H_ */
