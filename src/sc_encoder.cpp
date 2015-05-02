/*
 * sc_encoder.cpp
 *
 *  Created on: 2014/11/20
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cmath>
#include "sc_encoder.h"
#include "sc_utilities.h"
#include "bucket.h"

using namespace std;

namespace SC {
void SCEncoder::distribution(bool verbose) {
	size_t i, j;
	size_t base_u = 0, base_c = 0, base_p = 0;
	non_empty_bucket = 0;
	int hash, hash2;
	Bucket bk_tmp;
	size_t size1 = static_cast<size_t>(size);
	size_t size2 = (size1 * (size1 - 1)) >> 1;
	size_t kc = static_cast<size_t>(config.kc);
	int size3 = pow(kc,nc);
	size_t size4 = pow(kc,config.mc * nc);
	if(size4 >= INT_MAX) {
		cerr << "The size of this inverted index is too LARGE!" << endl;
		cerr << "kc=" << kc << ";mc=" << config.mc << ";nc=" << nc
				<< ";size4=" << size4 << endl;
		exit(EXIT_FAILURE);
	}
	size_t id;
	for(i = 0; i < size4; i++) {
		ivf.push_back(bk_tmp);
	}

	int uid[config.mc];
	for(i = 0; i < config.N; i++) {
		for(j = 0; j < config.mc; j++, base_u++) {
			uid[j] = 0;//config.kc * cid[base_u] + cid[base_u + config.mc];
			for(int k = 0; k < nc; k++) {
				uid[j] *= kc;
				uid[j] += cid[base_u + k * config.mc];
			}
		}

		// Set the bucket id
		id = vector_base(&uid[0],config.mc,size3);
		if(id != -1) {
			ivf[id].pid.insert(ivf[id].pid.end(),base_p++);
			ivf[id].codes.insert(ivf[id].codes.end(),
					&codes[base_c], &codes[base_c + config.mp]);
			ivf[id].L++;
			if(ivf[id].L == 1) non_empty_bucket++;
		}
		base_u += (config.mc);
		base_c += (config.mp);
	}

	if(verbose) {
		int sum = 0;
		int tmp = 0;
		double e = 0.0;
		double x;
		for(size_t i = 0; i < size4; i++) {
			tmp = ivf[i].L;
			if(tmp > 0) {
				//				cout << "ivf[" << i << "]:" << tmp << endl;
				x = 1.0 * config.N / tmp;
				e += log2(x) / x;
				sum += tmp;
			}
		}
		cout << sum << endl;
		cout << "The number of non-empty buckets: " << non_empty_bucket <<  endl;
		cout << "Entropy: " << e << endl;
		cout << "The maximum value of entropy: " << log2(non_empty_bucket)<< endl;
	}
}

void SCEncoder::output(
		const char * db_path,
		const char * db_prefix,
		bool verbose) {
	// Now we output the ivf structure to file
	char fname[256];
	sprintf(fname,"%s/%s_mr%d_ivf.edat_",db_path,db_prefix,nc);
	size_t size1 = static_cast<size_t>(size);
	size_t size2 = ((size1 * (size1 - 1)) >> 1);
	size_t kc = static_cast<size_t>(config.kc);
	int size3 = pow(kc,nc);
	size_t size4 = pow(size3,config.mc);
	if(size4 >= INT_MAX) {
		cerr << "The size of this inverted index is too LARGE!" << endl;
		cerr << "kc=" << kc << ";mc=" << config.mc << ";nc=" << nc
				<< ";size4=" << size4 << endl;
		exit(EXIT_FAILURE);
	}
	// Output the encoded data
#ifdef _WIN32
#else
	int fd = open(fname, O_RDWR | O_CREAT | O_TRUNC,(mode_t)0600); // file description
	if(fd < 0) {
		if(verbose)
			cerr << "Cannot open the file " << fname << endl;
		exit(1);
	}

	// The size of codes file
	size_t f_size = static_cast<size_t>(config.N)* static_cast<size_t>(config.mp)
															+ static_cast<size_t>(config.N + size4 + 2)
															* static_cast<size_t>(sizeof(int));
	int result = lseek(fd, f_size, SEEK_SET);
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

	unsigned char * fd_map = (unsigned char *)mmap(0, f_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (fd_map == MAP_FAILED) {
		close(fd);
		if(verbose)
			cerr << "Error mmapping the file" << endl;
		exit(1);
	}
	// Output the encoded data
	size_t bytes = 0;
	// The size of codes
	int N = config.N;
	// The first 8 bytes will be
	// the number of buckets
	memcpy(&fd_map[bytes],&non_empty_bucket,sizeof(int));
	bytes += sizeof(int);
	// and the size of the database
	memcpy(&fd_map[bytes],&N,sizeof(int));
	bytes += sizeof(int);

	// Output the buckets data
	size_t L, tmp, i, id, j, k;
	unsigned char c;
	vector<int> pid;
	vector<unsigned char> code2;
	for(i = 0; i < size4; i++) {
		L = ivf[i].L;
		// Write out the length of the bucket
		memcpy(&fd_map[bytes],&L,sizeof(int));
		bytes += sizeof(int);
		if(L > 0) {
			// Write the pid
			pid = ivf[i].pid;
			memcpy(&fd_map[bytes],&pid[0],L * sizeof(int));
			bytes += L *sizeof(int);

			// Write the code
			code2 = ivf[i].codes;
			memcpy(&fd_map[bytes],&code2[0],code2.size() * sizeof(unsigned char));
			bytes += code2.size() * sizeof(unsigned char);
			pid.clear();
			code2.clear();
		}
	}

	if (munmap(fd_map, f_size) == -1) {
		if(verbose)
			cerr << "Error un-mmapping the file" << endl;
	}
	close(fd);

	cout << "Read " << bytes << " byte(s) and wrote out " << f_size << " byte(s)" << endl;
#endif
}
} /* namespace PQLearn */
