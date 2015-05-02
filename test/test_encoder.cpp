/*
 * test_encoder.cpp
 *
 *  Created on: 2014/11/20
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include "encoder.h"
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
#include "sc_utilities.h"

using namespace std;
using namespace SC;

int param_k;

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
	static Encoder e;
};

// Global variables
Encoder EncoderTest::e;

TEST_F(EncoderTest, test1) {
	char filename1[256], filename2[256];
	sprintf(filename1,"./data/codebooks/cq_%d.ctr_",param_k);
	sprintf(filename2,"./data/codebooks/pq_%d.ctr_",param_k);
	e.load_codebooks(
			filename1,
			filename2,
			true);
	PQConfig config = e.get_config();
	EXPECT_EQ(1,config.mc);
	EXPECT_EQ(8,config.mp);
	EXPECT_EQ(128,config.dim);
	EXPECT_EQ(param_k,config.kc);
	EXPECT_EQ(256,config.kp);
}

TEST_F(EncoderTest, test2) {
	e.encode<float>("./data/sift/sift_base.fvecs",4,true);
	PQConfig config = e.get_config();
	EXPECT_EQ(1000000,config.N);
}

TEST_F(EncoderTest, test3) {
	e.distribution(true);
}

TEST_F(EncoderTest, test4) {
	char name[256];
	sprintf(name, "code_%d",param_k);
	e.output("./data/codebooks",name,true);
}

TEST_F(EncoderTest, DISABLED_test5) {
	e.set_mp(8);
	char filename[256];
	sprintf(filename, "./data/codebooks/code_%d_ivf.edat_",param_k);
	e.load_encoded_data(filename,true);
}

TEST_F(EncoderTest, DISABLED_test6) {
	e.rencode<unsigned char>("./data/sift/sift_base.bvecs",4,false);
}

TEST_F(EncoderTest, DISABLED_test7) {
	char name[256];
	sprintf(name, "code_%d_r",param_k);
	e.output2("./data/codebooks",name,true);
}

int main(int argc, char * argv[]) {
	/*
	 * The method is initializes the Google framework and must be called before RUN_ALL_TESTS
	 **/
	::testing::InitGoogleTest(&argc, argv);

	if(argc != 2) {
		cout << "Usage: " << argv[0] << " K" << endl;
		exit(EXIT_FAILURE);
	}

	param_k = atoi(argv[1]);

	/**
	 * RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	 * It's must be called only once in the code because multiple calls lead to conflicts and,
	 * therefore, are not supported.
	 */
	return RUN_ALL_TESTS();
}
