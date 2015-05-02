/*
 * pq_utilities.h
 *
 *  Created on: 2014/06/27
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef PQ_UTILITIES_H_
#define PQ_UTILITIES_H_

#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <iterator>
#include <vector>
#include <list>
#include <string>
#include <cstring>
#include <unordered_map>
#include <iterator>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cfloat>
#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
#include <utilities.h>

using namespace std;

namespace SC {
typedef unsigned short ushort;

/**
 * A projective function
 * f: [0,k) * [i,k) --> [0,k*(k-1) / 2)
 * @param k the base of the conversion
 * @param i
 * @param j
 */
inline int mapping(
		int k,
		int i,
		int j) {
	if(i == j || i >= k || j >= k || k < 0 || i < 0 || j < 0) return -1;
	if(i > j) return mapping(k,j,i);
	return i * k + j - (((i + 1) * (i + 2)) >> 1);
}

/**
 * A projective function
 * f^{-1}: [0,k*(k-1) / 2) --> [0,k) * [i,k)
 * @param n
 * @param k
 * @param i, j
 * @param tmp1, tmp2
 */
inline void resolve(
		int n,
		int k,
		int& i, int& j,
		float tmp1, float tmp2) {
	if(n < 0 || k < 0 || (n << 1) >= k * (k - 1)) {
		i = j = -1;
	}
	tmp1 = 1.0f * k - 0.5f;
	tmp2 = sqrt(tmp1 * tmp1 - 2.0f * n);
	i = static_cast<int>(tmp1 - tmp2);
	j = n - i * k + (((i + 1) * ( i + 2)) >> 1);
}

/**
 * A hash combiner from Boost C++ Library
 * http://www.boost.org/doc/libs/1_56_0/doc/html/hash/reference.html#boost.hash_combine
 * @param seed the seed to be computed
 * @param v the data to be hashed
 */
template<typename DataType>
inline void hash_combine(
		int& seed,
		DataType const& v) {
	seed ^= hash<DataType>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/**
 * An utility to get the file size
 * @param filename the name of the file
 * @return the size of the file as bytes
 */
inline size_t get_file_size(const char * filename) {
	if(filename == nullptr) exit(EXIT_FAILURE);
	size_t size, start;
	FILE * fp;
	fp = fopen(filename,"rb");
	if(fp == nullptr) exit(EXIT_FAILURE);
	start = ftell(fp);
	fseek(fp,0,SEEK_END);
	size = ftell(fp);
	fseek(fp,start,SEEK_SET);
	fclose(fp);
	return size;
}

/**
 * An utility to find the extension of file
 * @param filename the path to file
 * @param verbose enable verbose mode
 * @return the extension of file
 */
inline string file_extension(
		string filename,
		bool verbose) {
	string::size_type id = filename.rfind(".");
	if(id != string::npos) {
		return filename.substr(id+1);
	}
	return string("");
}

/**
 * Hash a vector of ushort components
 * @param code the vector
 * @param size the dimensionality of the vector
 * @return the hash value
 */
inline int vector_hash(
		ushort * code,
		int size) {
	int seed = 0;
	for(int i = 0; i < size; i++) {
		hash_combine<ushort>(seed,code[i]);
	}
	return seed;
}

/**
 * Convert a vector into a value of defined base
 * @param code the vector
 * @param size the dimensionality of the vector
 * @param base the base value
 * @return the conversion result as an integer
 */
inline int vector_base(
		ushort * code,
		int size,
		int base) {
	int seed = 0;
	for(int i = 0; i < size; i++) {
		seed *= base;
		seed += code[i];
	}
	return seed;
}

/**
 * Convert a vector into a value of defined base
 * @param code the vector
 * @param size the dimensionality of the vector
 * @param base the base value
 * @return the conversion result as a long integer
 */
inline size_t vector_base(
		int * code,
		int size,
		int base) {
	int seed = 0;
	for(int i = 0; i < size; i++) {
		seed *= base;
		seed += code[i];
	}
	return seed;
}

/**
 * Convert a number into a vector of defined base
 * @param v the input number
 * @param ut the vector
 * @param base the base value
 * @param size the dimensionality of the vector
 */
inline void vector_convert(
		int v,
		ushort *& ut,
		int base,
		int size) {
	if(base <= 0)
		exit(EXIT_FAILURE);
	int i = size - 1;
	while(i >= 0) {
		ut[i--] = v % base;
		v /= base;
	}
}

/**
 * Swap two integer by xor operations
 * @param a,b input number
 */
inline void swap_xor(int& a, int& b) {
	a ^= b;
	b ^= a;
	a ^= b;
}

/**
 * An utility to insert elements into a priority queue.
 * The priority queue is implemented by an array of float values.
 * @param heap the priority queue
 * @param hid the corresponding identifiers
 * @param data the value to be inserted
 * @param id the identifier of the data
 * @param pos the position of the data that is inserted
 * @param max_size the capacity of the heap
 * @verbose to see some log messages
 */
inline void pq_insert(
		float * heap,
		int * hid,
		float data,
		int id,
		int& pos,
		int max_size,
		bool verbose) {
	if(pos >= max_size || pos < 0) {
		if(verbose)
			cerr << "Heap's full or pos is negative" << endl;
		return;
	}

	heap[pos] = data;
	hid[pos] = id;
	pos++;
	int i = pos - 1, j = i >> 1;
	while(i > 0) {
		if(heap[i] < heap[j]) {
			// swapping data
			data = heap[j];
			heap[j] = heap[i];
			heap[i] = data;
			// swapping id
			swap_xor(hid[i],hid[j]);
		} else break;
		i = j;
		j = i >> 1;
	}
}

/**
 * Pop a value and its identifier from the priority queue
 * @param heap the priority queue
 * @param hid the identifiers
 * @param curpos the current capacity of the heap
 * @param verbose to enable verbose mode
 * @return the identifier of popped data
 */
inline int pq_pop(
		float * heap,
		int * hid,
		int& curpos,
		bool verbose) {
	if(curpos <= 0)
		exit(EXIT_FAILURE);
	int tmp = hid[0];
	float tmp2;
	heap[0] = heap[curpos - 1];
	hid[0] = hid[--curpos];
	int i = 0, j = i << 1, max;
	while(j + 1 < curpos) {
		if(heap[j] < heap[j+1]) {
			max = j;
		} else {
			max = j + 1;
		}

		if(heap[max] < heap[i]) {
			// swapping id
			swap_xor(hid[i],hid[max]);
			// swapping data
			tmp2 = heap[i];
			heap[i] = heap[max];
			heap[max] = tmp2;
			i = max;
		} else break;
		j = i << 1;
	}
	return tmp;
}

/**
 * An utility to insert elements into a priority queue.
 * The priority queue is implemented by an array of float values.
 * Each data in the priority queue has many identifiers.
 * @param heap the priority queue
 * @param hid the corresponding identifiers
 * @param hs the number of identifiers that assigned into each data
 * @param data the value to be inserted
 * @param id the identifier of the data
 * @param pos the position of the data that is inserted
 * @param max_size the capacity of the heap
 * @verbose to see some log messages
 */
inline void pq_insert_m(
		float * heap,
		int ** hid,
		int * id,
		int hs,
		float data,
		int &pos,
		int max_size,
		bool verbose) {
	if(pos >= max_size || pos < 0 || hs < 0) {
		cerr << "Heap's full or pos is negative:pos=" << pos << endl;
		return;
	}

	heap[pos] = data;
	for(int l = 0; l < hs; l++) {
		hid[l][pos] = id[l];
	}
	pos++;
	int i = pos - 1, j = i >> 1;
	while(i > 0) {
		if(heap[i] < heap[j]) {
			// swapping data
			data = heap[j];
			heap[j] = heap[i];
			heap[i] = data;
			// swapping id
			for(int l = 0; l < hs; l++) {
				swap_xor(hid[l][i],hid[l][j]);
			}
		} else break;
		i = j;
		j = i >> 1;
	}
}

/**
 * Pop a value and its identifier from the priority queue
 * @param heap the priority queue
 * @param hid the identifiers
 * @param hs the number of identifiers that are assigned into each data
 * @param curpos the current capacity of the heap
 * @param id the identifiers that are popped
 * @param v output value
 * @param verbose to enable verbose mode
 * @return the identifier of popped data
 */
inline void pq_pop_m(
		float * heap,
		int ** hid,
		int hs,
		int& curpos,
		int *& id,
		float& v,
		bool verbose) {
	if(curpos <= 0)
		return;

	for(int l = 0; l < hs; l++) {
		id[l] = hid[l][0];
	}
	v = heap[0];

	float tmp;
	heap[0] = heap[curpos - 1];
	for(int l = 0; l < hs; l++) {
		hid[l][0] = hid[l][curpos - 1];
	}
	--curpos;
	int i = 0, j = i << 1, max;
	while(j + 1 < curpos) {
		if(heap[j] < heap[j+1]) {
			max = j;
		} else {
			max = j + 1;
		}

		if(heap[max] < heap[i]) {
			// swapping id
			for(int l = 0; l < hs; l++) {
				swap_xor(hid[l][i],hid[l][max]);
			}
			// swapping data
			tmp = heap[i];
			heap[i] = heap[max];
			heap[max] = tmp;
			i = max;
		} else break;
		j = i << 1;
	}
}

inline void pq_insert2(float * heap, int * hid, int * hid2,
		float data, int id, int id2,
		int& pos, int max_size, bool verbose) {
	if(pos >= max_size || pos < 0) {
		if(verbose)
			cerr << "Heap's full or pos is negative:pos=" << pos << endl;
		return;
	}

	heap[pos] = data;
	hid[pos] = id;
	hid2[pos] = id2;
	pos++;
	int i = pos - 1, j = i >> 1;
	while(i > 0) {
		if(heap[i] < heap[j]) {
			// swapping data
			data = heap[j];
			heap[j] = heap[i];
			heap[i] = data;
			// swapping id
			swap_xor(hid[i],hid[j]);
			// swapping id2
			swap_xor(hid2[i],hid2[j]);
		} else break;
		i = j;
		j = i >> 1;
	}
}

inline void pq_pop2(float * heap, int * hid, int * hid2,
		int& curpos, int& id, int& id2, bool verbose) {
	if(curpos <= 0)
		return;
	id = hid[0];
	id2 = hid2[0];
	float tmp3;
	heap[0] = heap[curpos - 1];
	hid[0] = hid[curpos - 1];
	hid2[0] = hid2[--curpos];
	int i = 0, j = i << 1, max;
	while(j + 1 < curpos) {
		if(heap[j] < heap[j+1]) {
			max = j;
		} else {
			max = j + 1;
		}

		if(heap[max] < heap[i]) {
			// swapping id
			swap_xor(hid[i],hid[max]);
			swap_xor(hid2[i],hid2[max]);
			// swapping data
			tmp3 = heap[i];
			heap[i] = heap[max];
			heap[max] = tmp3;
			i = max;
		} else break;
		j = i << 1;
	}
}

inline void pq_insert4(
		float * heap,
		int * hid1, int * hid2, int * hid3, int * hid4,
		float data,
		int h1, int h2, int h3, int h4,
		int& pos, int max_size, bool verbose) {
	if(pos >= max_size || pos < 0) {
		if(verbose)
			cerr << "Heap's full or pos is negative:pos=" << pos << endl;
		return;
	}

	heap[pos] = data;
	hid1[pos] = h1;
	hid2[pos] = h2;
	hid3[pos] = h3;
	hid4[pos] = h4;
	pos++;
	int i = pos - 1, j = i >> 1;
	while(i > 0) {
		if(heap[i] < heap[j]) {
			// swapping data
			data = heap[j];
			heap[j] = heap[i];
			heap[i] = data;
			// swapping id
			swap_xor(hid1[i],hid1[j]);
			swap_xor(hid2[i],hid2[j]);
			swap_xor(hid3[i],hid3[j]);
			swap_xor(hid4[i],hid4[j]);
		} else break;
		i = j;
		j = i >> 1;
	}
}

inline void pq_insert6(
		float * heap,
		int * hid1, int * hid2, int * hid3, int * hid4,
		int * hid5, int * hid6,
		float data,
		int h1, int h2, int h3, int h4, int h5, int h6,
		int& pos, int max_size, bool verbose) {
	if(pos >= max_size || pos < 0) {
		if(verbose)
			cerr << "Heap's full or pos is negative:pos=" << pos << endl;
		return;
	}

	heap[pos] = data;
	hid1[pos] = h1;
	hid2[pos] = h2;
	hid3[pos] = h3;
	hid4[pos] = h4;
	hid5[pos] = h5;
	hid6[pos] = h6;
	pos++;
	int i = pos - 1, j = i >> 1;
	while(i > 0) {
		if(heap[i] < heap[j]) {
			// swapping data
			data = heap[j];
			heap[j] = heap[i];
			heap[i] = data;
			// swapping id
			swap_xor(hid1[i],hid1[j]);
			swap_xor(hid2[i],hid2[j]);
			swap_xor(hid3[i],hid3[j]);
			swap_xor(hid4[i],hid4[j]);
			swap_xor(hid5[i],hid5[j]);
			swap_xor(hid6[i],hid6[j]);
		} else break;
		i = j;
		j = i >> 1;
	}
}

inline void pq_pop4(
		float * heap,
		int * hid1, int * hid2, int * hid3, int * hid4,
		int& h1, int& h2, int &h3, int& h4,
		int& curpos, bool verbose) {
	if(curpos <= 0)
		return;
	h1 = hid1[0];
	h2 = hid2[0];
	h3 = hid3[0];
	h4 = hid4[0];
	float tmp;
	heap[0] = heap[curpos - 1];
	hid1[0] = hid1[curpos - 1];
	hid2[0] = hid2[curpos - 1];
	hid3[0] = hid3[curpos - 1];
	hid4[0] = hid4[--curpos];
	int i = 0, j = i << 1, max;
	while(j + 1 < curpos) {
		if(heap[j] < heap[j+1]) {
			max = j;
		} else {
			max = j + 1;
		}

		if(heap[max] < heap[i]) {
			// swapping id
			swap_xor(hid1[i],hid1[max]);
			swap_xor(hid2[i],hid2[max]);
			swap_xor(hid3[i],hid3[max]);
			swap_xor(hid4[i],hid4[max]);
			// swapping data
			tmp = heap[i];
			heap[i] = heap[max];
			heap[max] = tmp;
			i = max;
		} else break;
		j = i << 1;
	}
}

inline void pq_pop4_f(
		float * heap,
		int * hid1, int * hid2, int * hid3, int * hid4,
		int& h1, int& h2, int &h3, int& h4, float& v,
		int& curpos, bool verbose) {
	if(curpos <= 0)
		return;
	h1 = hid1[0];
	h2 = hid2[0];
	h3 = hid3[0];
	h4 = hid4[0];
	v = heap[0];
	float tmp;
	heap[0] = heap[curpos - 1];
	hid1[0] = hid1[curpos - 1];
	hid2[0] = hid2[curpos - 1];
	hid3[0] = hid3[curpos - 1];
	hid4[0] = hid4[--curpos];
	int i = 0, j = i << 1, max;
	while(j + 1 < curpos) {
		if(heap[j] < heap[j+1]) {
			max = j;
		} else {
			max = j + 1;
		}

		if(heap[max] < heap[i]) {
			// swapping id
			swap_xor(hid1[i],hid1[max]);
			swap_xor(hid2[i],hid2[max]);
			swap_xor(hid3[i],hid3[max]);
			swap_xor(hid4[i],hid4[max]);
			// swapping data
			tmp = heap[i];
			heap[i] = heap[max];
			heap[max] = tmp;
			i = max;
		} else break;
		j = i << 1;
	}
}

inline void pq_pop6(
		float * heap,
		int * hid1, int * hid2, int * hid3, int * hid4,
		int * hid5, int * hid6,
		int& h1, int& h2, int &h3, int& h4, int& h5, int& h6,
		int& curpos, bool verbose) {
	if(curpos <= 0)
		return;
	h1 = hid1[0];
	h2 = hid2[0];
	h3 = hid3[0];
	h4 = hid4[0];
	h5 = hid5[0];
	h6 = hid6[0];
	float tmp;
	heap[0] = heap[curpos - 1];
	hid1[0] = hid1[curpos - 1];
	hid2[0] = hid2[curpos - 1];
	hid3[0] = hid3[curpos - 1];
	hid4[0] = hid4[curpos - 1];
	hid5[0] = hid5[curpos - 1];
	hid6[0] = hid6[--curpos];
	int i = 0, j = i << 1, max;
	while(j + 1 < curpos) {
		if(heap[j] < heap[j+1]) {
			max = j;
		} else {
			max = j + 1;
		}

		if(heap[max] < heap[i]) {
			// swapping id
			swap_xor(hid1[i],hid1[max]);
			swap_xor(hid2[i],hid2[max]);
			swap_xor(hid3[i],hid3[max]);
			swap_xor(hid4[i],hid4[max]);
			swap_xor(hid5[i],hid5[max]);
			swap_xor(hid6[i],hid6[max]);
			// swapping data
			tmp = heap[i];
			heap[i] = heap[max];
			heap[max] = tmp;
			i = max;
		} else break;
		j = i << 1;
	}
}

/**
 * The linear search
 * @param data the input data
 * @param query the query data
 * @param N the size of the data
 * @param d the dimensionality of the vectors
 * @param verbose to enable the verbose mode
 * @return the identifier of the result
 */
template<typename DataType>
inline int linear_search(
		DataType * data,
		DataType * query,
		int N,
		int d,
		bool verbose) {
	float d_tmp, d_min = FLT_MAX;
	DataType * tmp = data;
	int i, ans = -1;
	for(i = 0; i < N; i++) {
		d_tmp = SimpleCluster::distance_l2_square<DataType>(query,tmp,d);
		if(d_min > d_tmp) {
			d_min = d_tmp;
			ans = i;
		}
		tmp += d;
	}
	if(verbose)
		cout << "Found at " << ans << " with distance = " << d_min << endl;
	return ans;
}

/**
 * The exact k-nn search
 * @param data the input data
 * @param query the query data
 * @param v_tmp a temporary array
 * @param heap a priority queue
 * @param hid the identifiers of the data in the priority queue
 * @param ans the result
 * @param N the size of the data
 * @param d the dimensionality of the vectors
 * @param k the number of nn to be retrieved
 * @param verbose to enable the verbose mode
 */
template<typename DataType>
inline void linear_knn(
		DataType * data,
		DataType * query,
		float * v_tmp,
		float * heap,
		int * hid,
		int * ans,
		int N,
		int d,
		int k,
		bool verbose) {
	DataType * tmp = data;
	int i, pos = 0;
	for(i = 0; i < N; i++) {
		v_tmp[i] = v_tmp[i+N] =
				SimpleCluster::distance_l2_square<DataType>(query,tmp,d);
		tmp += d;
	}
	nth_element(v_tmp+N,v_tmp+N+k-1,v_tmp+(N<<1));
	float d_tmp = v_tmp[N+k-1];
	for(i = 0; i < N; i++) {
		if(v_tmp[i] <= d_tmp) {
			pq_insert(heap,hid,v_tmp[i],i,pos,k,verbose);
		}
	}
	pos = k;
	for(i = 0; i < k; i++) {
		ans[i] = pq_pop(heap,hid,pos,verbose);
	}
}

/**
 * Strip matrix
 * Maybe using OpenBLAS will help the performance
 * @param data the data to be stripped
 * @param ans the stripped data
 * @param N the size of the input data
 * @param d the dimensionality of the input data
 * @param st,ed the interval to be stripped
 * @param verbose to enable the verbose mode
 */
template<typename DataType>
inline void strip_matrix(
		DataType * data,
		DataType *& ans,
		int N,
		int d,
		int st, int ed,
		bool verbose) {
	if(verbose)
		cout << "Strip components from " << st << " to " << ed << " of " << N << " vectors" << endl;
	if(st > ed || st < 0 || ed > d) {
		cerr << "Wrong parameters:st=" << st << ";ed=" << ed << ";d=" << d << endl;
		exit(EXIT_FAILURE);
	}
	size_t bs = ed - st, i, col, base = 0;

	SimpleCluster::init_array<DataType>(ans,bs * static_cast<size_t>(N));
	for(i = 0; i < static_cast<size_t>(N) * static_cast<size_t>(d); i++) {
		col = i % d;
		if(col < ed && col >= st) {
			ans[base++] = data[i];
		}
	}
	if(verbose)
		cout << "Stripped " << base << " components" << endl;
}
/**
 * Configuration structure
 * @param N the number of data
 * @param kc, kp the number of coarse quantizer and product quantizer's centers
 * @param mc, mp the number of partitions
 * @param dim the number of dimensions
 * @param w the number of bucket will be assigned into query when searching is performed
 * @param T the list's length
 * @param db_path path to the DB files
 * @param db_prefix the prefix of DB files
 */
typedef struct {
	int N, kc, kp, mc, mp, w, dim, T, L;
	char db_prefix[256];
	char db_path[256];
}PQConfig;

/**
 * Load data from file
 * @param filename the path to file
 * @param data the data
 * @param header the number of bytes of the header
 * @param d the number of dimensions
 * @param verbose enable verbose mode
 * @return the number of records that are loaded
 */
template<typename DataType>
inline int load_data(
		const char * filename,
		DataType *& data,
		int header,
		int d,
		bool verbose) {
#ifdef _WIN32
#else
	// Precheck
	if(d <= 0 || filename == nullptr) {
		if(verbose)
			cerr << "Defective data!" << endl;
		exit(1);
	}
	int fd = open(filename, O_RDONLY);
	if(fd < 0) {
		if(verbose)
			cerr << "Cannot open the file" << endl;
		exit(1);
	}

	size_t size = get_file_size(filename);
	unsigned char * mapped; // Read file as bytes
	/* Mapping the file */
	mapped = (unsigned char *)mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	if(mapped == MAP_FAILED) {
		if(verbose)
			cerr << "Cannot map the file" << endl;
		exit(1);
	}

	size_t i, base = header, count = 0;
	DataType f[1];
	size_t d1 = d * sizeof(DataType) + header;
	unsigned char uc, buf[sizeof(DataType)];
	size_t total_row = size / d1;
	data = (DataType *)::operator new(total_row * d * sizeof(DataType));
	if(verbose)
		cout << "We will load " << total_row << " vector(s)" << endl;

	/* Load data */
	while(count < total_row * d && base < size) {
		for(i = 0; i < sizeof(DataType); i++) {
			buf[i] = mapped[base++];
		}
		memcpy(f,buf,sizeof(DataType));
		data[count++] = f[0];
		if(count % d == 0) base += header;
	}

	if(verbose)
		cout << "The number of vectors: " << count / d << endl;

	if(munmap(mapped,size)  != 0) {
		cerr << "Cannot munmap file" << endl;
		exit(EXIT_FAILURE);
	}
	close(fd);

	if(verbose) {
		cout << "Last row:" << endl;
		for(i = (total_row - 1) * d; i < total_row  * d; i++) {
			cout << data[i] << " ";
		}
		cout << endl;
	}

	return count / d;
#endif
}

/**
 * Load a codebook from a binary file
 * @param filename the location of codebook files
 * @param config the configuration
 * @param codebook the codebook
 * @param type type=0 means codebook is a coarse codebook, otherwise it's a product codebook
 * @param verbose enable verbose mode
 */
template<typename DataType>
inline void load_codebook(
		const char * filename,
		PQConfig& config,
		DataType *& codebook,
		int type,
		bool verbose) {
	int bs;
	int k, m, row, col;
#ifdef _WIN32
#else
	int fd = open(filename, O_RDONLY);
	if(fd < 0) {
		if(verbose)
			cerr << "Cannot open the file" << endl;
		exit(1);
	}

	struct stat s;
	int status = fstat(fd, &s);
	if(status < 0) {
		if(verbose)
			cerr << "Cannot get statistics of file" << endl;
		exit(1);
	}

	size_t size = s.st_size;
	int count = size / sizeof(DataType);

	DataType * mapped;
	/* Mapping the file */
	mapped = (DataType *)mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	if(mapped == MAP_FAILED) {
		if(verbose)
			cerr << "Cannot map the file" << endl;
		exit(1);
	}

	config.dim = static_cast<int>(mapped[2]);
	if(type == 0) {
		config.kc = static_cast<int>(mapped[0]);
		k = config.kc;
		config.mc = static_cast<int>(mapped[1]);
		m = config.mc;
	} else if(type == 1) {
		config.kp = static_cast<int>(mapped[0]);
		k = config.kp;
		config.mp = static_cast<int>(mapped[1]);
		m = config.mp;
	} else {
		config.kc = static_cast<int>(mapped[0]);
		k = config.kc;
		config.mc = static_cast<int>(mapped[1]);
		m = config.mc;
		config.L = static_cast<int>(mapped[3]);
	}
	if(m == 0) {
		cerr << "Wrong data" << endl;
		exit(EXIT_FAILURE);
	}

	bs = config.dim / m;

	if(type == 0 || type == 1) {
		SimpleCluster::init_array<DataType>(codebook, k * config.dim);
		memcpy(codebook,mapped + 3,k * config.dim * sizeof(DataType));
	} else if(type == 3){
		SimpleCluster::init_array<DataType>(codebook, static_cast<size_t>(k) * (config.L + 1) * config.dim);
		memcpy(codebook,mapped + 4,static_cast<size_t>(k) * (config.L + 1) * config.dim * sizeof(DataType));
	} else if(type == 2){
		SimpleCluster::init_array<DataType>(codebook, (config.L + k) * config.dim);
		memcpy(codebook,mapped + 4,(config.L + k) * config.dim * sizeof(DataType));
	}

	if(munmap(mapped,size)  != 0) {
		cerr << "Cannot munmap file" << endl;
		exit(EXIT_FAILURE);
	}
	close(fd);
#endif

	if(verbose)
		cout << "Loaded a codebook at " << filename
		<< " with " << m << " sub codes" << endl;
}

/**
 * Save data into file and return the number of bytes that are written.
 * The input data in type DataType1 will be saved in DataType2.
 * @param filename output file
 * @param input input data array
 * @param row the number of rows
 * @param col the number of column
 * @param type type=true to use offset or type=false to disable offset
 * @param verbose true to enable verbose mode
 * @return Return the number of bytes that are written.
 */
template<typename DataType1, typename DataType2>
inline size_t save_data(
		const char * filename,
		DataType1 * input,
		size_t row,
		size_t col,
		bool type,
		bool verbose) {
	size_t size = 0;
#ifdef _WIN32
#else
	int fd = open(filename, O_RDWR | O_CREAT | O_TRUNC,(mode_t)0600); // file description
	if(fd < 0) {
		if(verbose)
			cerr << "Cannot open the file " << filename << endl;
		exit(1);
	}

	// The number of bytes to be written out.
	size = row * col * sizeof(DataType2);
	if(type)
		size += sizeof(int) * row;

	int result = lseek(fd, size, SEEK_SET);
	if (result == -1) {
		close(fd);
		if(verbose)
			cerr << "Error calling lseek() to 'stretch' the file" <<  endl;
		exit(1);
	}

	int status = write(fd, "", 1);
	if(status != 1) {
		if(verbose)
			cerr << "Cannot write to file" << endl;
		exit(1);
	}

	unsigned char * fd_map = (unsigned char *)mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (fd_map == MAP_FAILED) {
		close(fd);
		if(verbose)
			cerr << "Error mmapping the file" << endl;
		exit(1);
	}
	// Output the encoded data
	size_t bytes = 0;
	DataType2 tmp;

	size_t base = 0;
	for(int i = 0; i < row; i++) {
		if(type) {
			memcpy(&fd_map[bytes],&i,sizeof(int));
			bytes += sizeof(i);
		}
		for(int j = 0; j < col; j++) {
			tmp = static_cast<DataType2>(input[base++]);
			memcpy(&fd_map[bytes],&tmp,sizeof(DataType2));
			bytes += sizeof(tmp);
		}
	}

	if(verbose)
		cout << "Read " << bytes << " byte(s) from memory and saved "
		<< size << " byte(s) to disk" << endl;

	if (munmap(fd_map, size) == -1) {
		if(verbose)
			cerr << "Error un-mmapping the file" << endl;
	}
	close(fd);
#endif
	return size;
}

/**
 * Load data from file and convert into new data type
 * @param filename the path to file
 * @param data the data
 * @param header the number of bytes of the header
 * @param d the number of dimensions
 * @param verbose enable verbose mode
 * @return the number of records that are read
 */
template<typename DataType1, typename DataType2>
inline int load_and_convert_data(
		const char * filename,
		DataType2 *& data,
		int header,
		int d,
		bool verbose) {
#ifdef _WIN32
#else
	// Precheck
	if(d <= 0 || filename == nullptr) {
		if(verbose)
			cerr << "Defective data!" << endl;
		exit(1);
	}
	int fd = open(filename, O_RDONLY);
	if(fd < 0) {
		if(verbose)
			cerr << "Cannot open the file " << filename << endl;
		exit(1);
	}
	size_t size = get_file_size(filename);
	if(verbose)
		cout << "File size is " << size << endl;
	unsigned char * mapped; // Read file as bytes
	/* Mapping the file */
	mapped = (unsigned char *)mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	if(mapped == MAP_FAILED) {
		if(verbose)
			cerr << "Cannot map the file" << endl;
		exit(1);
	}

	size_t i, row = 0, base = header;
	size_t count = 0;
	DataType1 f[1];
	DataType2 f2;
	size_t d1 = d * sizeof(DataType1) + header;
	unsigned char buf[sizeof(DataType1)];
	size_t total_row = size / d1;
	data = (DataType2 *)::operator new(total_row * d * sizeof(DataType2));
	if(verbose)
		cout << "We will load " << total_row << " vector(s)" << endl;

	/* Load data */
	while(count < total_row * d && base < size) {
		for(i = 0; i < sizeof(DataType1); i++) {
			buf[i] = mapped[base++];
		}
		memcpy(f,buf,sizeof(DataType1));
		f2 = static_cast<DataType2>(f[0]);
		data[count++] = f2;
		if(count % d == 0) base += header;
	}

	if(verbose)
		cout << "The number of vectors: " << count / d << endl;

	if(munmap(mapped,size)  != 0) {
		cerr << "Cannot munmap file" << endl;
		exit(EXIT_FAILURE);
	}
	close(fd);

	return count / d;
#endif
}

/**
 * An utility to get the overview information of a files
 * @param data the data of the file
 * @param m,n the dimensionalities of the data
 */
template <typename DataType>
inline void get_statistic_of_data(
		DataType * data,
		int m, int n) {
	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
	//	omp_lock_t lock;
#endif
	double max = DBL_MIN, min = DBL_MAX, e1 = 0.0, e2 = 0.0, tmp;
	size_t p = static_cast<size_t>(m) / max_threads;
#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for private(tmp)
#endif
		for(int i0 = 0; i0 < max_threads; i0++) {
			size_t start = p * i0;
			size_t end = start + p;
			if(end > m || i0 ==  max_threads - 1) end = m;
			size_t base = start * n;
			for(int i = start; i < end; i++) {
				for(int j = 0; j < n; j++) {
					tmp = static_cast<double>(data[base++]);
					if(min > tmp) {
						min = tmp;
					}
					if(max < tmp) max = tmp;
					e1 += tmp;
					e2 += tmp * tmp;
				}
			}
		}
#ifdef _OPENMP
	}
#endif

	e1 /= (static_cast<size_t>(m) * static_cast<size_t>(n));

	cout << "Minimum component: " << min << endl;
	cout << "Maximum component: " << max << endl;
	cout << "Mean average: " << e1 << endl;
	cout << "Distribution: " << e2 / (static_cast<size_t>(m) * static_cast<size_t>(n)) - e1 * e1 << endl;
}

/**
 * Create a random data from a high dimensionality data
 * @param data the input data
 * @param subset the result data
 * @param N the size of the input data
 * @param d the dimensionality of the input data
 * @param k choose randomly one component from k components
 * @param verbose to enable verbose mode
 */
template<typename DataType>
inline void create_subset(
		DataType * data,
		DataType *& subset,
		size_t N,
		size_t d,
		int k,
		bool verbose) {
	if(k <= 1) return;
	size_t i, j, M = N / k, bs = d * k, r = N - M * k;
	DataType * tmp = data;
	SimpleCluster::init_array(subset,M * d);
	DataType * tmp2 = subset;
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dist(0,k - 1);
	for(i = 0; i < M - 1; i++) {
		j = dist(gen);
		memcpy(tmp2,tmp + j * d,d * sizeof(DataType));
		tmp2 += d;
		tmp += bs;
	}
	if(r == 0) {
		j = dist(gen);
	} else {
		uniform_int_distribution<int> dist2(0,r - 1);
		j = dist2(gen);
	}
	memcpy(tmp2,tmp + j * d,d * sizeof(DataType));
}
}

#endif /* PQ_UTILITIES_H_ */
