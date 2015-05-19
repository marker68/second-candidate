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

char config_file[256];

/**
 * Customized test case for testing
 */
class EncoderTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		ifstream input;
		input.open(config_file, ios::in);
		input >> cq_in >> mrq_in >> base_data >> codebook_out;
		input.close();
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
	static char cq_in[256], mrq_in[256], base_data[256], codebook_out[256];
};

// Global variables
SCEncoder EncoderTest::e(3);

char EncoderTest::cq_in[256];
char EncoderTest::mrq_in[256];
char EncoderTest::base_data[256];
char EncoderTest::codebook_out[256];

TEST_F(EncoderTest, test1) {
	e.load_codebooks(
			cq_in,
			mrq_in,
			true);
}

TEST_F(EncoderTest, test2) {
	e.encode<unsigned char>(base_data,4,false);
}

TEST_F(EncoderTest, test3) {
	e.distribution(true);
}

TEST_F(EncoderTest, test4) {
	e.output("./data/codebooks",codebook_out,true);
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
