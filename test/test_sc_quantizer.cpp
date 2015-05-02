/*
 * test_sc_quantizer.cpp
 *
 *  Created on: 2014/12/27
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cstring>
#include <gtest/gtest.h>
#include "sc_quantizer.h"
#include "sc_utilities.h"

using namespace std;
using namespace SC;

/**
 * Customized test case for testing
 */
class QuantizerTest : public ::testing::Test {
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
		::delete mrq;
		mrq = nullptr;
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static SCQuantizer<float> * mrq;
};

SCQuantizer<float> * QuantizerTest::mrq;

/**
 * Create product quantizer(m=8)
 */
TEST_F(QuantizerTest, test1) {
	mrq = new SCQuantizer<float>(128,8,256,4,1,true);
	mrq->load_data("./data/sift/sift_learn.fvecs",true);
}

TEST_F(QuantizerTest, test2) {
	char * filename1[] = { "./data/codebooks/sift_test1_cq.ctr_"};
	mrq->load_codebooks(filename1,true);
}

TEST_F(QuantizerTest, test3) {
	mrq->calc_residual_vector(true);
}

TEST_F(QuantizerTest, test4) {
	mrq->create_sub_quantizers(false);
}

TEST_F(QuantizerTest, test5) {
	cout << "Distortion is " << mrq->distortion(false) << endl;
}

TEST_F(QuantizerTest, test6) {
	mrq->output("./data/codebooks/sift_test1_mrq",true);
}

int main(int argc, char * argv[]) {
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);

	/*RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	It's must be called only once in the code because multiple calls lead to conflicts and,
	therefore, are not supported.
	*/
	return RUN_ALL_TESTS();
}
