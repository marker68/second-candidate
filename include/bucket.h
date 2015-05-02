/*
 * bucket.h
 *
 *  Created on: 2014/11/26
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef BUCKET_H_
#define BUCKET_H_

#include<iostream>
#include <vector>

using namespace std;
namespace SC {

/**
 * a struct to represent a bucket or a cell in the inverted file
 */
typedef struct {
	int L = 0; // the number of data in this cell
	vector<int> pid; // the identifiers of data in this cell
	vector<unsigned char> codes; // the codes of data in this cell
} Bucket;

typedef struct {
	int L = 0;
	vector<int> pid;
	vector<unsigned char> codes;
	int h, l;
} TreeBucket;

} /* namespace PQLearn */

#endif /* BUCKET_H_ */
