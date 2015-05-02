/*
 * sc_algorithm.h
 *
 *  Created on: 2014/12/10
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef SC_ALGORITHM_H_
#define SC_ALGORITHM_H_
#include <iostream>
#include <vector>
#include <algorithm>
#include <cblas.h>
#include "sc_utilities.h"

using namespace std;

namespace SC {
/**
 * Sorting with index tracking using STL's tools
 * @param st,ed pointers to specify the range of array to be sorted
 * @param id the id of the array
 */
inline void sort_id(float * st, float * ed, int * id) {
	size_t n = ed - st;
	vector<pair<int,float>> pv;

	float * tmp;
	int * it;
	for(tmp = st, it = id; tmp != ed; tmp++, it++) {
		pv.push_back(pair<int,float>(*it,*tmp));
	}

	sort(pv.begin(), pv.end(),
			[](const pair<int,float>& a, const pair<int,float>& b) -> bool
			{ return a.second < b.second ; });
	tmp = st;
	it = id;
	for(size_t i = 0; i < n; i++) {
		*tmp = pv[i].second;
		*it = pv[i].first;
		tmp++;
		it++;
	}
	pv.clear();
}

/**
 * nth_element with id tracking
 * @param st,ed pointers to specify the range of array to be sorted
 * @param id the id of the array
 */
inline void nth_element_id(float * st, float * ed, int * id, int r) {
	size_t n = ed - st;
	vector<pair<int,float>> pv;

	float * tmp;
	int * it;
	for(tmp = st, it = id; tmp != ed; tmp++, it++) {
		pv.push_back(pair<int,float>(*it,*tmp));
	}

	nth_element(pv.begin(), pv.begin() + r, pv.end(),
			[](const pair<int,float>& a, const pair<int,float>& b) -> bool
			{ return a.second < b.second ; });
	tmp = st;
	it = id;
	for(size_t i = 0; i < n; i++) {
		*tmp = pv[i].second;
		*it = pv[i].first;
		tmp++;
		it++;
	}
	pv.clear();
}

/**
 * Selection Sorting with index tracking
 * @param st,ed pointers to specify the range of array to be sorted
 * @param id the id of the array
 */
inline void select_sort_id(float * st, float * ed, int * id) {
	int N = ed - st;
	if(N < 1) return;
	int i, j = 0, pos, it;
	float tmp, min;
	for(i = 0; i < N - 1; i++) {
		min = st[i];
		pos = i;
		for(j = i + 1; j < N; j++) {
			if(st[j] <= min) {
				min = st[j];
				pos = j;
			}
		}
		if(pos > i) {
			it = id[pos];
			id[pos] = id[i];
			id[i] = it;
			tmp = st[pos];
			st[pos] = st[i];
			st[i] = tmp;
		}
	}
}

/**
 * Partitioning with index tracking
 * @param data the data to be partitioned
 * @param id the id of the array
 * @param pivot the pivot to partition array
 * @param N the size of the array data
 * @return the index of the separator
 */
inline int partition_id(float * data, float pivot, int N, int * id) {
	int i = 0, j = N - 1, it;
	float tmp;
	while(i <= j) {
		while(data[i] < pivot)
			i++;
		while(data[j] >= pivot)
			j--;
		if(i <= j) {
			tmp = data[i];
			data[i] = data[j];
			data[j] = tmp;
			it = id[i];
			id[i] = id[j];
			id[j] = it;
			i++;
			j--;
		}
	}

	return i;
}

/**
 * Quick Sorting with index tracking
 * @param st,ed pointers to specify the range of array to be sorted
 * @param id the id of the array
 */
inline void quick_sort_id(float * st, float * ed, int * id)  {
	int N = ed - st;
	if(N <= 5) {
		select_sort_id(st,ed,id);
		return;
	}
	// First, choose an appropriate pivot
	float pivot = st[N >> 1];
	// Choose the pivot as the median of left, right and (left+right)/2 elements
	if((pivot <= st[0]&& st[0] <= st[N-1])
			|| ((pivot >= st[0]&& st[0] >= st[N-1])))
		pivot = st[0];
	if((pivot <= st[N-1]&& st[N-1] <= st[0])
			|| ((pivot >= st[N-1]&& st[N-1] >= st[0])))
		pivot = st[N-1];

	// Then partition the array into two parts base on the pivot value.
	// The left part will contains members that are less than pivot,
	// the right part contains member that are greater than or equal pivot.
	int p = partition_id(st,pivot,N,id);
	if(p <= 0) {
		select_sort_id(st,ed,id);
		return;
	}

	// Finally, do the QuickSort recursively
	quick_sort_id(st,st + p,id);
	quick_sort_id(st + p,ed,id + p);
}

/**
 * Binary search
 */
inline void binary_search(float * data, float pivot,
		int& id, bool& found,int N,
		int st, int ed) {
	id = -1;
	found = false;
	if(N <= 0) {
		cerr << "Nothing to do" << endl;
		return;
	}
	if(st < 0) st = 0;
	if(ed >= N) ed = N - 1;

	int md;
	float tmp;
	if(ed - st > 1) {
		md = (st + ed) >> 1;
		tmp = data[md];
		if(tmp <= pivot) {
			binary_search(data,pivot,id,found,N,md,ed);
		} else {
			binary_search(data,pivot,id,found,N,st,md);
		}
	} else {
		id = st;
		if(data[st] == pivot) {
			found = true;
		} else {
			found = false;
		}
	}
}

/**
 * Calculate L2 distance between two vector with OpenBLAS
 */
inline float distance_l2(float * x, float * y, int n) {
	cblas_saxpy(n,-1,x,1,y,1);
	return cblas_snrm2(n,x,1);
}
}



#endif /* SC_ALGORITHM_H_ */
