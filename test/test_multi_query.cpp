/*
 * test_multi_query.cpp
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
#include <multi_query.h>

using namespace std;
using namespace SC;

char setting[256];

inline void create_log_dir(const char * log_path, time_t timer) {
	struct stat sb;
	char dir[256];
	if(stat(log_path, &sb) != 0) {
#ifdef _WIN32
		if(CreateDirectory(log_path, nullptr) == 0) {
#else
		if(mkdir(log_path,00777) == -1) {
#endif
			cerr << "Failed to create data folder search" << log_path << endl;
#ifdef _WIN32
			cerr << "Error code: " << GetLastError() << endl;
#endif
			exit(1);
		}
	}

	sprintf(dir,"%s/%d",log_path,static_cast<int>(timer));
#ifdef _WIN32
	if(CreateDirectory(dir, nullptr) == 0) {
#else
	if(mkdir(dir,00777) == -1) {
#endif
		cerr << "Failed to create data folder" << dir << endl;
#ifdef _WIN32
		cerr << "Error code: " << GetLastError() << endl;
#endif
		exit(1);
	}
}

/**
 * Customized test case for testing
 */
class QueryTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		strcpy(setting_path, setting);
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
	setting_path[256],
	base_path[256],
	query_path[256],
	cq_path[256],
	pq_path[256],
	code_path[256],
	log_path[256];
	static MultiQuery * worker;
	static float * data;
	static int d, w, T, N, kc, M;
};

// Global variables
char QueryTest::base_path[256];
char QueryTest::setting_path[256];
char QueryTest::query_path[256];
char QueryTest::cq_path[256];
char QueryTest::pq_path[256];
char QueryTest::code_path[256];
char QueryTest::log_path[256];
MultiQuery * QueryTest::worker;
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
	input.open(setting_path, ios::in);
	input >> base_path >> query_path
	>> cq_path >> pq_path >> code_path >> log_path
	>> d >> w >> T >> kc >> M;
	input.close();
}

TEST_F(QueryTest, test1) {
	worker = new MultiQuery();
}

TEST_F(QueryTest, test2) {
	worker->load_codebooks(
			cq_path,
			pq_path,true);
	worker->load_encoded_data(code_path,false);
	cout << "Entropy: " << worker->entropy(worker->get_full_size()) << endl;
	cout << "Max entropy:" << log2(worker->get_size()) << endl;
}

TEST_F(QueryTest, test3) {
	N = load_and_convert_data<unsigned char,float>(query_path,data,4,d,true);
//	worker->load_data<unsigned char>(base_path,4,true);
	worker->pre_compute1();
}

TEST_F(QueryTest, test4) {
	clock_t st, ed;
	time_t timer = time(NULL);
	create_log_dir(log_path,timer);
	double t = 0.0;
	int R, r[] = {1,10,100,1000,10000,100000};
	ofstream output;
	char filename[256];
	int * result, * it, * hid1, * hid2, * hid3, * hid4, * s1, * s2, * prebuck, * cache;
	float * v_tmp, * dist, * q, * tmp;
	bool * traversed;
	SimpleCluster::init_array(v_tmp,(kc << 1) + w);
	SimpleCluster::init_array(q, 2);
	SimpleCluster::init_array(it,kc<<1);
	SimpleCluster::init_array(hid1,w);
	SimpleCluster::init_array(hid2,w);
	SimpleCluster::init_array(hid3,w);
	SimpleCluster::init_array(hid4,w);
	SimpleCluster::init_array(s1,w);
	SimpleCluster::init_array(s2,w);
	SimpleCluster::init_array(prebuck,w);
	SimpleCluster::init_array(cache,w);
	SimpleCluster::init_array(traversed,kc * kc);
	double t1, t2;

	int sum = 0, e = 0;
	for(int i = 0; i < 3; i++) {
		t = 0.0;
		t1 = t2 = 0.0;
		R = r[i];
		tmp = data;
		sprintf(filename,"%s/%d/search_result_%d.txt",log_path,static_cast<int>(timer),R);
		output.open(filename,ios::out);
		for(int j = 0; j < M; j++) {
			st = clock();
			worker->search_multi2(
					tmp,
					v_tmp,dist,q,
					it,result,
					hid1,hid2,hid3,hid4,
					s1,s2,
					prebuck,cache,traversed,
					sum,R,w,T,kc,e,t1,t2,false,false);
			ed = clock();
			t += static_cast<double>(ed - st);
			EXPECT_FALSE(traversed[0]);
			for(int i = 0; i < (R>sum?sum:R); i++)
				output << result[i] << " ";
			if(sum < R)
				for(int i = 0; i < R - sum; i++)
					output << "-1 ";
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
		cout << "Precomputed in " <<
				1000.0f * t1 / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Sorted data in " <<
				1000.0f * t2 / CLOCKS_PER_SEC << "[ms]" << endl;
	}
}

TEST_F(QueryTest, DISABLED_test5) {
	clock_t st, ed;
	time_t timer = time(NULL);
	create_log_dir(log_path,timer);
	double t = 0.0;
	int R, r[] = {1,10,100,1000,10000,100000};
	ofstream output;
	char filename[256];
	int * result, * it, * hid1, * hid2, * hid3, * hid4, * s1, * s2, * prebuck, * cache;
	float * v_tmp, * dist, * q, * tmp;
	bool * traversed;
	SimpleCluster::init_array(v_tmp,(kc << 1) + w);
	SimpleCluster::init_array(q, 2);
	SimpleCluster::init_array(it,kc<<1);
	SimpleCluster::init_array(hid1,w);
	SimpleCluster::init_array(hid2,w);
	SimpleCluster::init_array(hid3,w);
	SimpleCluster::init_array(hid4,w);
	SimpleCluster::init_array(s1,w);
	SimpleCluster::init_array(s2,w);
	SimpleCluster::init_array(prebuck,w);
	SimpleCluster::init_array(cache,w);
	SimpleCluster::init_array(traversed,kc * kc);
	double t1, t2;

	int sum = 0, e = 0;
	for(int i = 0; i < 3; i++) {
		t = 0.0;
		t1 = t2 = 0.0;
		R = r[i];
		tmp = data;
		sprintf(filename,"%s/%d/search_result_%d.txt",log_path,static_cast<int>(timer),R);
		output.open(filename,ios::out);
		for(int j = 0; j < N; j++) {
			st = clock();
			worker->search_multi2(
					tmp,
					v_tmp,dist,q,
					it,result,
					hid1,hid2,hid3,hid4,
					s1,s2,
					prebuck,cache,traversed,
					sum,R,w,T,kc,e,t1,t2,true,false);
			ed = clock();
			t += static_cast<double>(ed - st);
			EXPECT_FALSE(traversed[0]);
			for(int i = 0; i < (R>sum?sum:R); i++)
				output << result[i] << " ";
			if(sum < R)
				for(int i = 0; i < R - sum; i++)
					output << "-1 ";
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
		cout << "Precomputed in " <<
				1000.0f * t1 / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Sorted data in " <<
				1000.0f * t2 / CLOCKS_PER_SEC << "[ms]" << endl;
	}
}

int main(int argc, char * argv[])
{
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);

	if(argc != 2) {
		cerr << "ERROR: Not enough arguments" << endl;
		cerr << "Usage: " << argv[0] << " /path/to/the/config/file" << endl;
		exit(EXIT_FAILURE);
	}
	// Parse the arguments
	sprintf(setting,"%s",argv[1]);

	/**
	 * RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	 * It's must be called only once in the code because multiple calls lead to conflicts and,
	 * therefore, are not supported.
	 **/
	return RUN_ALL_TESTS();
}
