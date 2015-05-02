/*
 * sc_query.cpp
 *
 *  Created on: 2014/12/19
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <algorithm>
#include "sc_query.h"

using namespace std;

namespace SC
{

SCQuery::SCQuery() : PQQuery::PQQuery()
{
	nc = 2;
}

SCQuery::SCQuery(int _nc) : PQQuery::PQQuery()
{
	nc = _nc;
}

SCQuery::~SCQuery()
{
}

/**
 * Load the encoded data from binary file
 * @param filename path to the encoded data file
 * @param verbose enable verbose mode
 */
int SCQuery::load_encoded_data(const char * filename, bool verbose) {
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

	size_t f_size = s.st_size; // The size of file

	unsigned char * mapped, * temp;
	/* Mapping the file */
	mapped = (unsigned char *)mmap(0, f_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if(mapped == MAP_FAILED) {
		if(verbose)
			cerr << "Cannot map the file " << filename << endl;
		exit(1);
	}

	size_t l, i, j, k, count = 0, tmp;
	size_t base_pid = 0, base_code = 0;
	size_t size1 =  static_cast<size_t>(size);
	size_t size2 = pow(size1,nc);
	temp = mapped;

	// Read the number of buckets and the size of database
	memcpy(&not_empty,temp,sizeof(int));
	cout << "The number of non empty buckets: " << not_empty << endl;
	temp += sizeof(int);
	memcpy(&(config.N),temp,sizeof(int));
	temp += sizeof(int);

	// Memory allocation
	// SimpleCluster::init_array<int>(bits,size >> 5 + 1); // The list of non-empty buckets as a bit set
	SimpleCluster::init_array(L,size2); // Only restore the length of non-empty buckets
	SimpleCluster::init_array(pid,config.N);
	codes = (unsigned char *)::operator new(static_cast<size_t>(config.N) *
			static_cast<size_t>(config.mp) * sizeof(unsigned char));

	for(i = 0; i < size2; i++) {
		// Read the length of the bucket
		memcpy(&l,temp,sizeof(int));
		temp += sizeof(int);
		//		if(l > 0) {
		//			tmp = i & 32;
		//			bits[i>>5] |= (1 << (31 - tmp));
		//			L[count++] = l;
		//		}
		if(i > 0) L[i] = l + L[i-1];
		else L[i] = l;
		if(verbose)
			cout << "This bucket contains " << l << " vector(s)" << endl;

		if(l > 0) {
			// Read the pid
			memcpy(pid + base_pid,temp,l * sizeof(int));
			temp += l * sizeof(int);
			base_pid += l;

			// Read the codes
			memcpy(codes + base_code,temp,l * config.mp * sizeof(unsigned char));
			temp += l * config.mp * sizeof(unsigned char);
			base_code += l * config.mp;
		}
	}


	if(munmap(mapped,f_size) != 0) {
		cerr << "Cannot munmap file data" << endl;
		exit(EXIT_FAILURE);
	}
	close(fd);

	cout << "Read " << config.N << " data  from " << filename << endl;

	return config.N;
#endif
}

} /* namespace PQLearn */
