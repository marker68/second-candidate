/*
 * evaluation.h
 *
 *  Created on: 2014/06/30
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef EVALUATION_H_
#define EVALUATION_H_

#include <iostream>
#include <fstream>
#include <utilities.h>
#include "evaluation.h"
#include "sc_utilities.h"

using namespace std;

namespace SC {
class Evaluation {
public:
	inline void load_groundtruth(const char *,int *&, int,int&,bool);
	inline void calc_recall(int *, int *, int, int, int&, bool);
	inline void calc_comb_recall(int **, int *, int, int, int, int&, bool);
};

/**
 * Load the ground truth data
 * @param filename path to the ground truth file
 * @param d the number of components each row
 * @param N the number of queries
 * @param verbose enable verbose mode
 */
inline void Evaluation::load_groundtruth(
		const char * filename,
		int *& gt,
		int d,
		int& N,
		bool verbose) {
	int * tmp;
	N = load_data(filename,tmp,4,d,verbose);
	SimpleCluster::init_array<int>(gt,N);
	for(int i = 0; i < N; i++) {
		gt[i] = *tmp;
		tmp += d;
	}
}

/**
 * Calculate the recall@R
 * @param result search result
 * @param groundtruth the ground truth data
 * @param R
 * @param N the number of queries
 * @param recall the recall@R
 * @param verbose enable verbose mode
 */
inline void Evaluation::calc_recall(
		int * result,
		int * groundtruth,
		int R,
		int N,
		int& recall,
		bool verbose) {
	recall = 0;
	int gt;
	int * tmp = result;
	for(int i = 0; i < N; i++) {
		//		tmp = result;
		gt = groundtruth[i];
		bool found = false;
		for(int j = 0; j < R; j++) {
			if(tmp[j] == gt) {
				found = true;
				break;
			}
		}
		if(found) recall++;
		if(/*found &&*/ verbose) {
			cout << recall << endl;
			cout << "GT: " << gt << endl;
			cout << "RESULT: ";
			for(int j = 0; j < R; j++) {
				cout << tmp[j] << " ";
			}
			cout << endl;
		}
		tmp += R;
	}
	if(verbose)
		cout << "Recall@" << R << " = " << recall << endl;
}

/**
 * Calculate the recall@R
 * @param result search result
 * @param groundtruth the ground truth data
 * @param R
 * @param N the number of queries
 * @param recall the recall@R
 * @param verbose enable verbose mode
 */
inline void Evaluation::calc_comb_recall(
		int ** results,
		int * groundtruth,
		int R,
		int N,
		int S,
		int& recall,
		bool verbose) {
	recall = 0;
	int gt;
	int ** tmp = (int **)::operator new(S * sizeof(int *));
	for(int i = 0; i < S; i++) {
		tmp[i] = results[i];
	}
	for(int i = 0; i < N; i++) {
		gt = groundtruth[i];
		bool found = false;
		for(int j = 0; j < S; j++) {
			for(int k = 0; k < R; k++) {
				if(tmp[j][k] == gt) {
					found = true;
					break;
				}
			}
		}
		if(found) recall++;
		if(verbose) {
			cout << recall << endl;
			cout << "GT: " << gt << endl;
			cout << endl;
		}
		for(int j = 0; j < S; j++) {
			tmp[j] += R;
		}
	}
	if(verbose)
		cout << "Recall@" << R << " = " << recall << endl;
}
}

#endif /* EVALUATION_H_ */
