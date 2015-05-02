/*
 * query.cpp
 *
 *  Created on: 2014/10/15
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cmath>
#include <query.h>

using namespace std;

namespace SC {

/**
 * The constructor
 */
PQQuery::PQQuery() {
	cq = nullptr;
	pq = nullptr;
	L = nullptr;
	pid = nullptr;
	codes =  nullptr;
	norm_c = nullptr;
	norm_r = nullptr;
	dot_cr = nullptr;
	diff_qr = nullptr;
	diff_qc = nullptr;
	real_dist = nullptr;
	raw_data = nullptr;
}

/**
 * The destructor
 */
PQQuery::~PQQuery() {
	::delete cq;
	::delete pq;
	::delete L;
	::delete pid;
	cq = nullptr;
	pq = nullptr;
	L = nullptr;
	pid = nullptr;
	codes =  nullptr;
	::delete norm_c;
	::delete norm_r;
	::delete dot_cr;
	::delete diff_qr;
	::delete diff_qc;
	::delete real_dist;
	::delete raw_data;
	norm_c = nullptr;
	norm_r = nullptr;
	dot_cr = nullptr;
	diff_qr = nullptr;
	diff_qc = nullptr;
	real_dist = nullptr;
	raw_data = nullptr;
}

/**
 * Load the codebooks from file
 */
void PQQuery::load_codebooks(
		const char * cq_path,
		const char *  pq_path,
		bool verbose) {
	load_codebook<float>(cq_path,config,cq,0,verbose);
	load_codebook<float>(pq_path,config,pq,1,verbose);
	size = static_cast<int>(pow(config.kc,config.mc));
	cout << "--> Settings: (kc,mc,kp,mp)=" << config.kc << " "
			<< config.mc << " " << config.kp << " " << config.mp << endl;
}

/**
 * Load the encoded data from binary file
 * @param filename path to the encoded data file
 * @param verbose enable verbose mode
 */
int PQQuery::load_encoded_data(const char * filename, bool verbose) {
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

	int l;
	size_t i, j, k, count = 0, tmp;
	size_t base_pid = 0, base_code = 0;
	temp = mapped;

	// Read the number of buckets and the size of database
	memcpy(&not_empty,temp,sizeof(int));
	cout << "The number of non empty buckets: " << not_empty << "/" << size << endl;
	temp += sizeof(int);
	memcpy(&(config.N),temp,sizeof(int));
	temp += sizeof(int);

	// Memory allocation
	// SimpleCluster::init_array<int>(bits,size >> 5 + 1); // The list of non-empty buckets as a bit set
	SimpleCluster::init_array(L,size); // Only restore the length of non-empty buckets
	SimpleCluster::init_array(pid,config.N);
	cout << "We will load " << config.N << " elements" << endl;
	codes = (unsigned char *)::operator new(static_cast<size_t>(config.N) *
			static_cast<size_t>(config.mp) * sizeof(unsigned char));

	for(i = 0; i < size; i++) {
		// Read the length of the bucket
		memcpy(&l,temp,sizeof(int));
		temp += sizeof(int);
		if(verbose)
			cout << "This bucket contains " << l << " vector(s)" << endl;

		if(i > 0) L[i] = l + L[i-1];
		else L[i] = l;

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


double PQQuery::entropy(size_t size) {
	double e = 0.0;
	double N = config.N, l = L[0], x;
	x = N/l;
	if(L[0] > 0)
		e = log2(x) / x;
	for(int i = 1; i < size; i++) {
		l = L[i] - L[i-1];
		if(l > 0.0) {
			x = N / l;
			e += log2(x) / x;
		}
	}
	return e;
}

int PQQuery::get_size() {
	return not_empty;
}

int PQQuery::get_full_size() {
	return size;
}
} /* namespace PQLearn */
