/*
 * test_openblas.cpp
 *
 *  Created on: 2014/11/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <cblas.h>
#include <stdio.h>
#include <gtest/gtest.h>
#include <random>
#include <iostream>
#include <numeric>
#include <utilities.h>

using namespace std;

/**
 * Customized test case for testing
 */
class BLASTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		d1 = (float *)::operator new(64 * sizeof(float));
		d2 = (float *)::operator new(64 * sizeof(float));
		random_device rd;
		mt19937 gen(rd());
		normal_distribution<float> dis(0.0f,10000.0f);

		for(int i = 0; i < 64; i++) {
			d1[i] = dis(gen);
			d2[i] = dis(gen);
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
	static float * d1, * d2;
};

float * BLASTest::d1;
float * BLASTest::d2;

TEST_F(BLASTest, test1) {
	float d;
	for(int i = 0; i < 16384 * 2; i++) {
		d = inner_product(d1,d1+64,d2,0.0);
	}
}

TEST_F(BLASTest, test2) {
	float d;
	for(int i = 0; i < 16384 * 2; i++) {
		d = cblas_sdot(64,d1,1,d2,1);
	}
}

TEST_F(BLASTest, test3) {
	float v1,v2;
	v1 = inner_product(d1,d1+64,d2,0.0);
	v2 = cblas_sdot(64,d1,1,d2,1);
	EXPECT_LT(fabs(v1-v2),128.0);
}

TEST_F(BLASTest, test4) {
	float d;
	for(int i = 0; i < 16384 * 2; i++) {
		d = SimpleCluster::distance_l2_square(d1,d2,64);
	}
}

TEST_F(BLASTest, test5) {
	float d;
	for(int i = 0; i < 16384 * 2; i++) {
		cblas_saxpy(64,-1,d1,1,d2,1);
		d = cblas_snrm2(64,d2,1);
	}
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


