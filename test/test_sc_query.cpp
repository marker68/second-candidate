/*
 * test_sc_query.cpp
 *
 *  Created on: 2014/12/15
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
	db_prefix[256],
	cq_in[256],
	mrq_in[256],
	code_data[256];
	static SCQuery * worker;
	static float * data;
	static int d, w, T, N, kc, M;
};

// Global variables
char QueryTest::base_dir[256];
char QueryTest::query_path[256];
char QueryTest::db_prefix[256];
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

/**
 * Load configuration
 */
TEST_F(QueryTest, test0) {
	ifstream input;
	input.open(config_file, ios::in);
	input >> base_dir >> query_path >> db_prefix
	>> d >> w >> T >> kc >> M;
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
	N = load_and_convert_data<unsigned char,float>(query_path,data,4,d,true);
//	worker->load_data<unsigned char>("./data/sift/sift_base.bvecs",4,true);
	worker->pre_compute1();
}

TEST_F(QueryTest, test4) {
	clock_t st, ed;
	time_t timer = time(NULL);
	char data_folder[256];
	sprintf(data_folder,"%s/search",base_dir);
	struct stat sb;
if(stat(data_folder, &sb) != 0) {
#ifdef _WIN32
if(CreateDirectory(data_folder, nullptr) == 0) {
#else
	if(mkdir(data_folder,00777) == -1) {
#endif
		cerr << "Failed to create data folder search" << endl;
#ifdef _WIN32
		cerr << "Error code: " << GetLastError() << endl;
#endif
		exit(1);
	}
}
sprintf(data_folder,"%s/search/%d",base_dir,static_cast<int>(timer));
#ifdef _WIN32
if(CreateDirectory(data_folder, nullptr) == 0) {
#else
	if(mkdir(data_folder,00777) == -1) {
#endif
		cerr << "Failed to create data folder" << endl;
#ifdef _WIN32
		cerr << "Error code: " << GetLastError() << endl;
#endif
		exit(1);
	}
	double t = 0.0;
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
	for(int i = 0; i < 5; i++) {
		t = 0.0;
		R = r[i];
		tmp = data;
		sprintf(filename,"%s/search/%d/search_result_%d.txt",base_dir,static_cast<int>(timer),R);
		output.open(filename,ios::out);
		for(int j = 0; j < N; j++) {
			st = clock();
			worker->search_mr_ivf(
					tmp,
					v_tmp,dist,
					it,result,
					hid1,hid2,hid3,hid4,
					s1,s2,
					prebuck,cache,traversed,
					sum,R,w,T,M,false,false);
			ed = clock();
			t += static_cast<double>(ed - st);
			EXPECT_TRUE(traversed[0]);
			EXPECT_TRUE(!traversed[1]);
			for(int i = 0; i < (R>sum?sum:R); i++) {
				output << result[i] << " ";
			}
			if(sum < R)
				for(int i = 0; i < R - sum; i++) {
					output << "-1 ";
				}
			output << endl;
			tmp += d;
			::delete result;
			result = nullptr;
			::delete dist;
			dist =  nullptr;
		}

		output.close();
		cout << "Finished search@" << R << " in " <<
				1000.0f * t / CLOCKS_PER_SEC << "[ms]" << endl;
	}
}

TEST_F(QueryTest, DISABLED_test5) {
	clock_t st, ed;
	time_t timer = time(NULL);
	char data_folder[256];
	sprintf(data_folder,"%s/search",base_dir);
	struct stat sb;
if(stat(data_folder, &sb) != 0) {
#ifdef _WIN32
if(CreateDirectory(data_folder, nullptr) == 0) {
#else
	if(mkdir(data_folder,00777) == -1) {
#endif
		cerr << "Failed to create data folder search" << endl;
#ifdef _WIN32
		cerr << "Error code: " << GetLastError() << endl;
#endif
		exit(1);
	}
}
sprintf(data_folder,"%s/search/%d",base_dir,static_cast<int>(timer));
#ifdef _WIN32
if(CreateDirectory(data_folder, nullptr) == 0) {
#else
	if(mkdir(data_folder,00777) == -1) {
#endif
		cerr << "Failed to create data folder" << endl;
#ifdef _WIN32
		cerr << "Error code: " << GetLastError() << endl;
#endif
		exit(1);
	}
	double t = 0.0;
	int M = kc >> 4;
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
	t = 0.0;
	tmp = data;
	sprintf(filename,"%s/search/%d/search_result_%d.txt",base_dir,static_cast<int>(timer),1);
	output.open(filename,ios::out);
	for(int j = 0; j < N; j++) {
		st = clock();
		worker->search_mr_ivf(
				tmp,
				v_tmp,dist,
				it,result,
				hid1,hid2,hid3,hid4,
				s1,s2,
				prebuck,cache,traversed,
				sum,1,w,T,M,true,false);
		ed = clock();
		t += static_cast<double>(ed - st);
		for(int i = 0; i < 1; i++) {
			output << result[i] << " ";
		}
		output << endl;
		tmp += d;
		::delete result;
		result = nullptr;
		::delete dist;
		dist =  nullptr;
	}

	output.close();
	cout << "Finished search@1 in " <<
			1000.0f * t / CLOCKS_PER_SEC << "[ms]" << endl;
}

int main(int argc, char * argv[]) {
	/**
	 * The method is initializes the Google framework and must be called before RUN_ALL_TESTS
	 **/
	::testing::InitGoogleTest(&argc, argv);

	if(argc != 2) {
		cout << "Usage: " << argv[0] << " config_file" << endl;
		exit(EXIT_FAILURE);
	}

	snprintf(config_file, 256, "%s",argv[1]);

	/**
	 * RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	 * It's must be called only once in the code because multiple calls lead to conflicts and,
	 * therefore, are not supported.
	 * */
	return RUN_ALL_TESTS();
}
