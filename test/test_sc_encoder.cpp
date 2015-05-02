/*
 * test_sc_encoder.cpp
 *
 *  Created on: 2014/11/20
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include "sc_encoder.h"
#include <iostream>
#include <fstream>
#include <random>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <gtest/gtest.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif
#include "sc_encoder.h"
#include "sc_utilities.h"

using namespace std;
using namespace SC;

/**
 * Customized test case for testing
 */
class EncoderTest : public ::testing::Test {
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
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static SCEncoder e;
};

// Global variables
SCEncoder EncoderTest::e(3);

TEST_F(EncoderTest, test1) {
	e.load_codebooks(
			"./data/codebooks/sift_test64_8_cq.ctr_",
			"./data/codebooks/sift_test64_8_pq.ctr_",
			true);
	PQConfig config = e.get_config();
	EXPECT_EQ(1,config.mc);
	EXPECT_EQ(8,config.mp);
	EXPECT_EQ(128,config.dim);
	EXPECT_EQ(64,config.kc);
	EXPECT_EQ(256,config.kp);
}

TEST_F(EncoderTest, test2) {
	e.encode<unsigned char>("./data/sift/sift_base.bvecs",4,false);
}

TEST_F(EncoderTest, test3) {
	e.distribution(true);
}

TEST_F(EncoderTest, test4) {
	e.output("./data/codebooks","sift_test64_8",true);
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
