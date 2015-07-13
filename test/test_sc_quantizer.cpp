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


char config_file[256];

/**
 * Customized test case for testing
 */
class QuantizerTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		ifstream input;
		input.open(config_file, ios::in);
		input >> learn_data >> cq_out >> mrq_out
		>> d >> m >> kc >> k >> type;
		input.close();
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
	static int d, m, kc, k, type;
	static char learn_data[256], cq_out[256], mrq_out[256];
};

SCQuantizer<float> * QuantizerTest::mrq;
int QuantizerTest::d;
int QuantizerTest::m;
int QuantizerTest::kc;
int QuantizerTest::k;
int QuantizerTest::type;
char QuantizerTest::learn_data[256];
char QuantizerTest::cq_out[256];
char QuantizerTest::mrq_out[256];

/**
 * Create product quantizer(m=8)
 */
TEST_F(QuantizerTest, test1) {
	mrq = new SCQuantizer<float>(d,m,kc,4,k,true);
	mrq->load_data(learn_data,true);
	mrq->set_type(type);
}

TEST_F(QuantizerTest, test2) {
	char * filename1[] = { cq_out};
	mrq->load_codebooks(filename1,true);
}

TEST_F(QuantizerTest, test3) {
	mrq->calc_residual_vector(true);
}

TEST_F(QuantizerTest, test4) {
	mrq->create_sub_quantizers(false);
}

TEST_F(QuantizerTest, test5) {
	if(mrq->get_type() == 2)
		mrq->distortion(false);
	cout << "Distortion is " << mrq->error << endl;
}

TEST_F(QuantizerTest, test6) {
	mrq->output(mrq_out,true);
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
