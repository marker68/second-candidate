/*
 * test_utilities.cpp
 *
 *  Created on: 2014/10/03
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cstring>
#include <gtest/gtest.h>
#include "sc_utilities.h"

using namespace std;
using namespace SC;

/**
 * Customized test case for testing
 */
class PQUtilTest : public ::testing::Test {
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
	static char filename[256];
	static float * data0;
	static int * i_data0;
	static unsigned char * data1;
};

char PQUtilTest::filename[256];
float * PQUtilTest::data0;
int * PQUtilTest::i_data0;
unsigned char * PQUtilTest::data1;

TEST_F(PQUtilTest, test1) {
	strcpy(filename,"./data/siftsmall/siftsmall_base.fvecs");
	EXPECT_EQ(10000, load_data<float>(filename,data0,4,128,false));
}

TEST_F(PQUtilTest, test2) {
	EXPECT_EQ(0,data0[0]);
	get_statistic_of_data<float>(data0,10000,128);
}

TEST_F(PQUtilTest, test3) {
	strcpy(filename,"./data/sift/sift_groundtruth.ivecs");
	EXPECT_EQ(10000, load_data<int>(filename,i_data0,4,100,false));
	get_statistic_of_data<int>(i_data0,10000,100);
}

TEST_F(PQUtilTest, test4) {
	EXPECT_EQ(string("bvecs"),file_extension(string("./data/test.bvecs"),true));
}

TEST_F(PQUtilTest, test5) {
	EXPECT_NE(string("bvecs"),file_extension(string("bvecs"),true));
}

TEST_F(PQUtilTest, test6) {
	size_t size = save_data<float,unsigned char>(
			"./data/siftsmall/siftsmall_base.bvecs",
			data0,10000,128,true,true);
	EXPECT_EQ(1320000, size);
}

TEST_F(PQUtilTest, test7) {
	strcpy(filename,"./data/siftsmall/siftsmall_base.bvecs");
	EXPECT_EQ(10000, load_data<unsigned char>(filename,data1,4,128,false));
}

TEST_F(PQUtilTest, test8) {
	get_statistic_of_data<unsigned char>(data1,10000,128);
}

TEST_F(PQUtilTest, test9) {
	float * ans;
	strip_matrix(data0,ans,10000,128,0,3,true);
	EXPECT_EQ(35,ans[2]);
	::delete ans;
	ans = nullptr;
}

TEST_F(PQUtilTest, test10) {
	float * query = data0 + 128 * 256;
	EXPECT_EQ(256,linear_search<float>(data0,query,10000,128,true));
}

TEST_F(PQUtilTest, test11) {
	float tmp1, tmp2;
	int i,j;
	resolve(0,256,i,j,tmp1,tmp2);
	EXPECT_EQ(0,i);
	EXPECT_EQ(1,j);
	resolve(254,256,i,j,tmp1,tmp2);
	EXPECT_EQ(0,i);
	EXPECT_EQ(255,j);
	resolve(255,256,i,j,tmp1,tmp2);
	EXPECT_EQ(1,i);
	EXPECT_EQ(2,j);
	resolve(32639,256,i,j,tmp1,tmp2);
	EXPECT_EQ(254,i);
	EXPECT_EQ(255,j);
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
