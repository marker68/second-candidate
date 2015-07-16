/*
 * query_sc.cpp
 *
 *  Created on: 2015/07/16
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstring>
#include <gtest/gtest.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif
#include <sc_query.h>

using namespace std;
using namespace SC;

char config_file[256];
int qid;

/**
 * Customized test case for testing
 */
class QueryTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
		::delete worker;
		worker = nullptr;
		::delete data;
		data = nullptr;
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static char
	base_dir[256],
	query_path[256],
	cq_in[256],
	mrq_in[256],
	code_data[256];
	static SCQuery * worker;
	static float * data;
	static int d, w, T, N, kc, M, offset;
};

// Global variables
char QueryTest::base_dir[256];
char QueryTest::query_path[256];
char QueryTest::cq_in[256];
char QueryTest::mrq_in[256];
char QueryTest::code_data[256];
SCQuery * QueryTest::worker;
float * QueryTest::data;
int QueryTest::d;
int QueryTest::w;
int QueryTest::T;
int QueryTest::N;
int QueryTest::kc;
int QueryTest::M;
int QueryTest::offset;

/**
 * Load configuration
 */
TEST_F(QueryTest, test0) {
	ifstream input;
	input.open(config_file, ios::in);
	input >> base_dir >> query_path >> cq_in >> mrq_in >> code_data
	>> d >> w >> T >> kc >> M >> offset;
	input.close();
}

TEST_F(QueryTest, test1) {
	worker = new SCQuery();
}

TEST_F(QueryTest, test2) {
	worker->load_codebooks(
			cq_in,
			mrq_in,
			true);
	worker->load_encoded_data(code_data,false);
	size_t size = worker->get_full_size();
	size_t size2 = size * size;
	cout << "Maximum non-empty cells: " << size2 - size << endl;
	cout << "Entropy: " << worker->entropy(size2) << endl;
	cout << "Max entropy:" << log2(worker->get_size()) << endl;
}

TEST_F(QueryTest, test3) {
	N = load_data<float>(query_path,data,offset,d,true);
	worker->pre_compute1();
}

TEST_F(QueryTest, test4) {
	struct stat sb;
	int R, r[] = {1,10,100,1000,10000,100000};
	ofstream output;
	char filename[256];
	int * result, * it, * hid1, * hid2, * hid3, * hid4, * s1, * s2, * prebuck, * cache;
	float * v_tmp, * dist, * tmp;
	bool * traversed;
	SimpleCluster::init_array(v_tmp,kc + w);
	SimpleCluster::init_array(it,kc);
	SimpleCluster::init_array(hid1,w);
	SimpleCluster::init_array(hid2,w);
	SimpleCluster::init_array(hid3,w);
	SimpleCluster::init_array(hid4,w);
	SimpleCluster::init_array(s1,w);
	SimpleCluster::init_array(s2,w);
	SimpleCluster::init_array(prebuck,w);
	SimpleCluster::init_array(cache,w);
	SimpleCluster::init_array(traversed,kc * kc);
	memset(traversed,0,kc * kc);
	traversed[0] = true;

	int sum = 0;
	for(int i = 0; i < 3; i++) {
		R = r[i];
		tmp = data + qid * d;
		sprintf(filename,"%s/sc_%d.txt",base_dir,R);
		output.open(filename,ios::out);
			worker->search_mr_ivf(
					tmp,
					v_tmp,dist,
					it,result,
					hid1,hid2,hid3,hid4,
					s1,s2,
					prebuck,cache,traversed,
					sum,R,w,T,M,false,false);
			EXPECT_TRUE(traversed[0]);
			EXPECT_TRUE(!traversed[1]);
			for(int i = 0; i < (R>sum?sum:R); i++) {
				output << result[i] << endl;
                cout << result[i] << endl;
			}
			if(sum < R)
				for(int i = 0; i < R - sum; i++) {
					output << "-1 " << endl;
				}
			::delete result;
			result = nullptr;
			::delete dist;
			dist =  nullptr;

		output.close();
	}
}

int main(int argc, char * argv[]) {
	/**
	 * The method is initializes the Google framework and must be called before RUN_ALL_TESTS
	 **/
	::testing::InitGoogleTest(&argc, argv);

	if(argc != 3) {
		cout << "Usage: " << argv[0] << " config_file qid" << endl;
		exit(EXIT_FAILURE);
	}

	snprintf(config_file, 256, "%s",argv[1]);
    qid = atoi(argv[2]);
	/**
	 * RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	 * It's must be called only once in the code because multiple calls lead to conflicts and,
	 * therefore, are not supported.
	 * */
	return RUN_ALL_TESTS();
}
