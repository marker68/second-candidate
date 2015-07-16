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

char config[256];

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
		input.open(config,ios::in);
		input >> base >> cq_path >> pq_path >> output_path >> output_name >> offset;
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
	static Encoder e;
	static char base[256], cq_path[256], pq_path[256], output_path[256], output_name[256];
	static int offset;
};

// Global variables
Encoder EncoderTest::e;
int EncoderTest::offset;
char EncoderTest::base[256];
char EncoderTest::cq_path[256];
char EncoderTest::pq_path[256];
char EncoderTest::output_path[256];
char EncoderTest::output_name[256];

TEST_F(EncoderTest, test1) {
	e.load_codebooks(
			cq_path,
			pq_path,
			true);
	PQConfig config = e.get_config();
}

TEST_F(EncoderTest, test2) {
	e.encode<unsigned char>(base,offset,true);
}

TEST_F(EncoderTest, test3) {
	e.distribution(true);
}

TEST_F(EncoderTest, test4) {
	e.output(output_path,output_name,true);
}

int main(int argc, char * argv[]) {
	/*
	 * The method is initializes the Google framework and must be called before RUN_ALL_TESTS
	 **/
	::testing::InitGoogleTest(&argc, argv);

	if(argc != 2) {
		cout << "Usage: " << argv[0] << " config_file" << endl;
		exit(EXIT_FAILURE);
	}

	snprintf(config,256,"%s",argv[1]);

	/**
	 * RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	 * It's must be called only once in the code because multiple calls lead to conflicts and,
	 * therefore, are not supported.
	 */
	return RUN_ALL_TESTS();
}
