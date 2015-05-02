/*
 * test_evaluation.cpp
 *
 *  Created on: 2014/10/09
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cstring>
#include <cstdio>
#include <climits>
#include <cerrno>
#include <gtest/gtest.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif
#include "sc_utilities.h"
#include "evaluation.h"

using namespace std;
using namespace SC;

/**
 * Global variables
 */
char path[256];

/**
 * A helper to get the list of sub-directories
 * @param dir_path the location of directory
 * @param list list of dir_path's sub-directories
 */
void get_list_dir(const char * dir_path, vector<string>& list) {
#ifdef _WIN32
#else
	DIR * dir;
	struct dirent * dirp;
	if((dir = opendir(dir_path)) == nullptr) {
		cerr << "Error(" << errno << ") while opening " << dir_path << endl;
		exit(1);
	}

	struct stat st;
	while((dirp = readdir(dir)) != nullptr) {
		lstat(dirp->d_name,&st);
		if(strcmp(dirp->d_name,".") != 0
				&& strcmp(dirp->d_name,"..") != 0) { // only if it's a directory
			list.push_back(string(dirp->d_name));
		}
	}
	closedir(dir);
#endif
	cout << "Got " << list.size() << " directories" << endl;
}

/**
 * Customized test case for testing
 */
class RecallTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		strcpy(setting_path, "./data/config/evaluation.txt");
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static char setting_path[256],
	base_dir[256],
	groundtruth_path[256];
	static int d;
	static int N;
	static int M;
	static Evaluation e;
	static int * gt;
};

// Global variables
char RecallTest::groundtruth_path[256];
char RecallTest::setting_path[256];
char RecallTest::base_dir[256];
int RecallTest::d;
Evaluation RecallTest::e;
int * RecallTest::gt;
int RecallTest::N;
int RecallTest::M;

/**
 * Load configuration
 */
TEST_F(RecallTest, test0) {
	ifstream input;
	input.open(setting_path, ios::in);
	input >> groundtruth_path
	>> d >> M;
	input.close();
}

TEST_F(RecallTest, test1) {
	e.load_groundtruth(groundtruth_path,gt,d,N,true);
}

TEST_F(RecallTest, test2) {
	ifstream input;
	vector<string> list_dir;
	get_list_dir(path,list_dir);
	int r[] = {1,10,100,1000,10000};
	int R;

	// Do the statistical stuff
	vector<string>::iterator it = list_dir.begin();
	vector<string>::iterator ie = list_dir.end();
	char gt_path[256];
	e.load_groundtruth(groundtruth_path,gt,d,N,true);
	while(it != ie) {
		cout << "Experiment was hold at " << *it << endl;
		for(int i = 0; i < 5; i++) {
			R = r[i];
			int * result;
			SimpleCluster::init_array(result,N*R);
			char result_path[256];
			sprintf(result_path, "%s/%s/search_result_%d.txt",path,it->c_str(),R);
			input.open(result_path,ios::in);
			int tmp;
			int count = 0;
			while(count < N * R) {
				input >> tmp;
				result[count++] = tmp;
			}
			input.close();
			int recall;
			e.calc_recall(result,gt,R,M,recall,false);
			cout << "Recall@" << R << " = " << 100.0f * recall / M << endl;
			::delete result;
			result = nullptr;
		}
		++it;
	}
	::delete gt;
	gt = nullptr;
}

int main(int argc, char * argv[])
{
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);
	if(argc != 2) {
		cerr << "Usage: " << argv[0] << " path_to_evaluation_folder" << endl;
		exit(EXIT_FAILURE);
	}
	sprintf(path,"%s",argv[1]);

	/*RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	It's must be called only once in the code because multiple calls lead to conflicts and,
	therefore, are not supported.
	*/
	return RUN_ALL_TESTS();
}


