/*
 * test_quantizer.cpp
 *
 *  Created on: 2014/09/04
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cstring>
#include <gtest/gtest.h>
#include "quantizer.h"
#include "sc_utilities.h"

using namespace std;
using namespace SC;

int param_k;

/**
 * Customized test case for testing
 */
class QuantizerTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		d = 128;
		m = 8;
		k = param_k;
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
		::delete cq;
		cq = nullptr;
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static char filename[256];
	static int d, m, k;
	static PQQuantizer<float> * cq;
};

char QuantizerTest::filename[256];
int QuantizerTest::d;
int QuantizerTest::m;
int QuantizerTest::k;
PQQuantizer<float> * QuantizerTest::cq;

TEST_F(QuantizerTest, test1) {
	strcpy(filename,"./data/sift/sift_base.fvecs");
	PQQuantizer<float> pq(d,m,k,4,true);
	pq.load_data(filename,true);
	EXPECT_EQ(1000000,pq.size());
}

/**
 * Create coarse quantizer(m=1)
 */
TEST_F(QuantizerTest, test2) {
	strcpy(filename,"./data/sift/sift_base.fvecs");
	cq = new PQQuantizer<float>(d,1,k,4,false);
	cq->load_data(filename,true);
	EXPECT_EQ(1000000,cq->size());
}

TEST_F(QuantizerTest, test3) {
	cq->create_sub_quantizers(false);
}

TEST_F(QuantizerTest, test4) {
	cout << "Distortion is " << cq->distortion(false) << endl;
}

TEST_F(QuantizerTest, test5) {
	sprintf(filename, "./data/codebooks/cq_%d",k);
	cq->output(filename,true);
}

/**
 * Create product quantizer(m=8)
 */
TEST_F(QuantizerTest, test6) {
	cq->set_params(1,m,256);
	EXPECT_EQ(1000000,cq->size());
}

TEST_F(QuantizerTest, test7) {
	sprintf(filename, "./data/codebooks/cq_%d.ctr_",k);
	char * filename1[] = {filename};
	cq->load_codebooks(filename1,true);
}

TEST_F(QuantizerTest, test8) {
	cq->calc_residual_vector(true);
}

TEST_F(QuantizerTest, test9) {
	cq->create_sub_quantizers(false);
}

TEST_F(QuantizerTest, test10) {
	cout << "Distortion is " << cq->distortion(false) << endl;
}

TEST_F(QuantizerTest, test11) {
	sprintf(filename, "./data/codebooks/pq_%d",k);
	cq->output(filename,true);
}

/**
 * Create product quantizer(m=8)
 */
TEST_F(QuantizerTest, DISABLED_test12) {
	cq->set_params(1,m,k);
	EXPECT_EQ(100000,cq->size());
}

TEST_F(QuantizerTest, DISABLED_test13) {
	sprintf(filename, "./data/codebooks/pq_%d",k);
	char * filename1[] = { filename};
	cq->load_codebooks(filename1,true);
}

TEST_F(QuantizerTest, DISABLED_test14) {
	cq->calc_residual_vector(true);
}

TEST_F(QuantizerTest, DISABLED_test15) {
	cq->create_sub_quantizers(false);
}

TEST_F(QuantizerTest, DISABLED_test16) {
	cout << "Distortion is " << cq->distortion(false) << endl;
}

TEST_F(QuantizerTest, DISABLED_test17) {
	sprintf(filename, "./data/codebooks/rq_%d",k);
	cq->output(filename,true);
}

int main(int argc, char * argv[]) {
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);

	if(argc != 2) {
		cout << "Usage: " << argv[0] << " K" << endl;
		exit(EXIT_FAILURE);
	}

	param_k = atoi(argv[1]);

	/*RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	It's must be called only once in the code because multiple calls lead to conflicts and,
	therefore, are not supported.
	*/
	return RUN_ALL_TESTS();
}
