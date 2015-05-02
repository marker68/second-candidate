/*
 * encoder.cpp
 *
 *  Created on: 2014/11/20
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cmath>
#include "encoder.h"
#include "sc_utilities.h"
#include "bucket.h"

using namespace std;

namespace SC {

Encoder::Encoder() {
	cq = nullptr;
	pq = nullptr;
	cid =  nullptr;
	codes = nullptr;
	old_mp = 0;
	L = nullptr;
	pid = nullptr;
	L = nullptr;
	pid = nullptr;
}

Encoder::~Encoder() {
	::delete cq;
	::delete pq;
	::delete cid;
	::delete codes;
	::delete L;
	::delete pid;
	cq = nullptr;
	pq = nullptr;
	cid =  nullptr;
	codes = nullptr;
	ivf.clear();
}

void Encoder::load_codebooks(
		const char * cq_path,
		const char * pq_path,
		bool verbose) {
	load_codebook<float>(cq_path,config,cq,0,verbose);
	load_codebook<float>(pq_path,config,pq,1,verbose);
	size = static_cast<int>(pow(config.kc,config.mc));
	if(verbose)
		cout << "The size of ivf:" << size << endl;
	cout << "--> Settings: (kc,mc,kp,mp)=" << config.kc << " "
			<< config.mc << " " << config.kp << " " << config.mp << endl;
}

void Encoder::load_encoded_data(
		const char * filename,
		bool verbose) {
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

	int i, j, k;
	size_t base_pid = 0, base_c = 0;
	temp = mapped;
	int not_empty, p, p1, l;

	// Read the number of buckets and the size of database
	memcpy(&non_empty_bucket,temp,sizeof(int));
	cout << "The number of non empty buckets: " << non_empty_bucket << endl;
	temp += sizeof(int);
	memcpy(&(config.N),temp,sizeof(int));
	temp += sizeof(int);

	// Memory allocation
	SimpleCluster::init_array(cid,static_cast<size_t>(config.N)
			* static_cast<size_t>(config.mc)); // Only restore the length of non-empty buckets
	cout << "Old size: " << old_mp << "[byte(s)]" << endl;
	SimpleCluster::init_array(L,size);
	SimpleCluster::init_array(pid,static_cast<size_t>(config.N));

	ushort c[config.mc];

	for(i = 0; i < size; i++) {
		// Read the length of the bucket
		memcpy(&l,temp,sizeof(int));
		L[i] = l;
		temp += sizeof(int);
		if(verbose)
			cout << "This bucket contains " << l << " vector(s)" << endl;
		p1 = i;

		for(j = config.mc - 1; j >= 0; --j) {
			c[j] = p1 % config.kc;
			p1 = (p1 - c[j]) / config.kc;
		}

		memcpy(&pid[base_pid],temp,l * sizeof(int));
		base_pid += l;
		temp += l * sizeof(int);

		for(j = 0; j < l; j++) {
			for(k = 0; k < config.mc; k++) {
				cid[base_c++] = c[k];
			}
		}
		temp += l * old_mp * sizeof(unsigned char);
		if(verbose)
			cout << "Read " << l << " vector(s)" << endl;
	}


	if(munmap(mapped,f_size) != 0) {
		cerr << "Cannot munmap file data" << endl;
		exit(EXIT_FAILURE);
	}
	close(fd);

	cout << "Read " << config.N << " data  from " << filename << endl;
#endif
}

PQConfig Encoder::get_config() {
	return config;
}

void Encoder::statistic(bool detail) {
	ofstream op;
	op.open("./ivf.txt",ios::out);
	int i, j, base = 0;
	for(i = 0; i < size; i++) {
		op << "Bucket " << i << ":" << ivf[i].L << endl;
		base = 0;
		for(j = 0; j < ivf[i].L; j++) {
			op << ivf[i].pid[j] << ":";
			for(int k = 0; k < config.mp; k++) {
				op << static_cast<int>(ivf[i].codes[base++]) << " ";
			}
			op << endl;
		}
	}
	op.close();
}

void Encoder::distribution(bool verbose) {
	size_t i, j;
	size_t base_u = 0, base_c = 0, base_p = 0;
	non_empty_bucket = 0;
	int hash;
	Bucket bk_tmp;
	ivf.clear();
	for(i = 0; i < size; i++) {
		ivf.push_back(bk_tmp);
	}
	for(i = 0; i < config.N; i++) {
		// Check whether if the bucket was created or not?
		hash = vector_base(&cid[base_u],config.mc,config.kc);
		// The bucket
		if(hash < size) {
			// Set the bucket id
			ivf[hash].pid.insert(ivf[hash].pid.end(),base_p++);
			ivf[hash].codes.insert(ivf[hash].codes.end(),
					&codes[base_c], &codes[base_c + config.mp]);
			ivf[hash].L++;
			if(ivf[hash].L == 1) non_empty_bucket++;
		}
		base_u += config.mc;
		base_c += config.mp;
	}

	if(verbose) {
		cout << "The number of non-empty buckets: " << non_empty_bucket <<  endl;
		int sum = 0;
		for(i = 0; i < size; i++) {
			cout << "ivf[" << i << "]:" << ivf[i].L << endl;
			sum += ivf[i].L;
		}
		cout << sum << endl;
	}
}

void Encoder::output(
		const char * db_path,
		const char * db_prefix,
		bool verbose) {
	// Now we output the ivf structure to file
	char fname[256];
	sprintf(fname,"%s/%s_ivf.edat_",db_path,db_prefix);
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
	size_t f_size = static_cast<size_t>(config.N)
					* static_cast<size_t>(config.mp)
					+ static_cast<size_t>(config.N + size + 2)
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
	int L, tmp, i, j;
	unsigned char c;
	vector<int> pid;
	vector<unsigned char> code2;
	for(i = 0; i < size; i++) {
		L = ivf[i].L;
		// Write out the length of the bucket
		memcpy(&fd_map[bytes],&L,sizeof(int));
		bytes += sizeof(int);

		// Write the pid
		pid = ivf[i].pid;
		for(j = 0; j < L; j++) {
			tmp = pid[j];
			memcpy(&fd_map[bytes],&tmp,sizeof(int));
			bytes += sizeof(int);
		}

		// Write the code
		code2 = ivf[i].codes;
		if(code2.size() != L * config.mp) {
			cerr << "Wrong data" << endl;
			exit(EXIT_FAILURE);
		}
		for(j = 0; j < code2.size(); j++) {
			c = code2[j];
			memcpy(&fd_map[bytes],&c,sizeof(unsigned char));
			bytes += sizeof(unsigned char);
		}
		pid.clear();
		code2.clear();
	}

	if (munmap(fd_map, f_size) == -1) {
		if(verbose)
			cerr << "Error un-mmapping the file" << endl;
	}
	close(fd);

	cout << "Read " << bytes << " byte(s) and wrote out " << f_size << " byte(s)" << endl;
#endif
}

void Encoder::output2(
		const char * db_path,
		const char * db_prefix,
		bool verbose) {
	// Now we output the ivf structure to file
	char fname[256];
	sprintf(fname,"%s/%s_ivf.edat_",db_path,db_prefix);
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
	size_t f_size = static_cast<size_t>(config.N)
											* static_cast<size_t>(config.mp)
											+ static_cast<size_t>(config.N + size + 2) * static_cast<size_t>(sizeof(int));
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
	int l, tmp, i, j;
	unsigned char c;
	int * _pid = pid;
	unsigned char * _code2 = codes;
	for(i = 0; i < size; i++) {
		l = L[i];
		// Write out the length of the bucket
		memcpy(&fd_map[bytes],&l,sizeof(int));
		bytes += sizeof(int);

		// Write the pid
		memcpy(&fd_map[bytes],_pid,l * sizeof(int));
		bytes += l * sizeof(int);
		_pid += l;

		// Write the code
		memcpy(&fd_map[bytes],_code2, l * config.mp * sizeof(unsigned char));
		bytes += l * config.mp * sizeof(unsigned char);
		_code2 += l * config.mp;
		if(verbose)
			cout << "Exported " << l << " vector(s)" << endl;
	}

	if (munmap(fd_map, f_size) == -1) {
		if(verbose)
			cerr << "Error un-mmapping the file" << endl;
	}
	close(fd);

	cout << "Read " << bytes << " byte(s) and wrote out " << f_size << " byte(s)" << endl;
#endif
}

void Encoder::set_mp(int old) {
	old_mp = old;
}

void Encoder::set_size(size_t _size) {
	size = _size;
}
} /* namespace PQLearn */
