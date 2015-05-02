/*
 * test_algorithm.cpp
 *
 *  Created on: 2014/10/03
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cstring>
#include <cmath>
#include <random>
#include <gtest/gtest.h>
#include <utilities.h>
#include "sc_algorithm.h"

using namespace std;
using namespace SC;

/**
 * Customized test case for testing
 */
class AlgorithmTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		N = 1<<20;
		SimpleCluster::init_array(data,N);
		SimpleCluster::init_array(id,N);
		SimpleCluster::init_array(data2,N);
		SimpleCluster::init_array(id2,N);
		SimpleCluster::init_array(data3,N);
		SimpleCluster::init_array(id3,N);
		// For generating random numbers
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<float> real_dis(0.0, 400000.0);

		for(size_t i = 0; i < N; i++) {
			data3[i] = data2[i] = data[i] = real_dis(gen);
			id3[i] = id2[i] = id[i] = i;
		}
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
	static float * data, * data2, * data3;
	static int * id, * id2, * id3;
	static size_t N;
};

float * AlgorithmTest::data;
float * AlgorithmTest::data2;
float * AlgorithmTest::data3;
size_t AlgorithmTest::N;
int * AlgorithmTest::id;
int * AlgorithmTest::id2;
int * AlgorithmTest::id3;
bool * td = (bool *)::operator new(1<<28);

TEST_F(AlgorithmTest, test1) {
	float arr[] = { 1.0, 3.2, 2.4, 2.3};
	int id[] = {0,1,2,3};
	sort_id(arr,arr+4,id);
	EXPECT_EQ(1.0,arr[0]);
	EXPECT_LT(fabs(2.3-arr[1]),1e-05);
	EXPECT_EQ(3,id[1]);
}

TEST_F(AlgorithmTest, test2) {
	sort_id(data,data+N,id);
}

TEST_F(AlgorithmTest, test3) {
	quick_sort_id(data2,data2+N,id2);
}

TEST_F(AlgorithmTest, test4) {
	for(int i = 0; i < N; i++) {
		EXPECT_EQ(data2[i],data[i]) << "data differ at " << i << endl;
		EXPECT_TRUE((id2[i] == id[i])
				|| (id2[i] != id[i] && (data[i] = data[i + 1] || data[i] == data[i-1]))) << "indices differ at " << i << endl;
	}
}

TEST_F(AlgorithmTest, test5) {
	memset(td,0,1<<28);
}

int main(int argc, char * argv[]) {
	/*
	 * The method is initializes the Google framework and must be called before RUN_ALL_TESTS
	 **/
	::testing::InitGoogleTest(&argc, argv);

	/**
	 * RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	 * It's must be called only once in the code because multiple calls lead to conflicts and,
	 * therefore, are not supported.
	 */
	return RUN_ALL_TESTS();
}
