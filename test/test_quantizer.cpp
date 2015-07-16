/*
 * quantizer_test.cpp
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

char config[256];

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
		input.open(config, ios::in);
		input >> filename >> cq_path >> pq_path >>
		d >> c >> m >> k >> offset >> type;
		input.close();
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
	static char cq_path[256], pq_path[256];
	static int d, c, m, k, offset, type;
	static PQQuantizer<unsigned char> * cq;
};

char QuantizerTest::filename[256];
char QuantizerTest::cq_path[256];
char QuantizerTest::pq_path[256];
int QuantizerTest::d;
int QuantizerTest::c;
int QuantizerTest::m;
int QuantizerTest::k;
int QuantizerTest::type;
int QuantizerTest::offset;
PQQuantizer<unsigned char> * QuantizerTest::cq;

/**
 * Create coarse quantizer(m=1)
 */
TEST_F(QuantizerTest, test1) {
	cq = new PQQuantizer<unsigned char>(d,c,k,offset,true);
	cq->set_type(type);
	cq->load_data(filename,true);
}

TEST_F(QuantizerTest, test2) {
	cq->create_sub_quantizers(false);
}

TEST_F(QuantizerTest, test3) {
	if(cq->get_type() == 2)
		cq->distortion(false);
	cout << "Distortion is " << cq->error << endl;
}

TEST_F(QuantizerTest, test4) {
	cq->output(cq_path,true);
}

/**
 * Create product quantizer(m=8)
 */
TEST_F(QuantizerTest, test5) {
	cq->set_params(1,m,256);
}

TEST_F(QuantizerTest, test6) {
	char tmp[256];
	snprintf(tmp,256,"%s.ctr_",cq_path);
	char * filename1[] = {tmp};
	cq->load_codebooks(filename1,true);
}

TEST_F(QuantizerTest, test7) {
	cq->calc_residual_vector(true);
}

TEST_F(QuantizerTest, test8) {
	cq->create_sub_quantizers(false);
}

TEST_F(QuantizerTest, test9) {
	if(cq->get_type() == 2)
		cq->distortion(false);
	cout << "Distortion is " << cq->error << endl;
}

TEST_F(QuantizerTest, test10) {
	cq->output(pq_path,true);
}

int main(int argc, char * argv[]) {
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);

	if(argc != 2) {
		cout << "Usage: " << argv[0] << " config_file" << endl;
		exit(EXIT_FAILURE);
	}

	snprintf(config,256,"%s",argv[1]);

	/*RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	It's must be called only once in the code because multiple calls lead to conflicts and,
	therefore, are not supported.
	*/
	return RUN_ALL_TESTS();
}
