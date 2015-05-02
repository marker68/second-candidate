/*
 * test_hash.cpp
 *
 *  Created on: 2014/10/03
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cstring>
#include <random>
#include <gtest/gtest.h>
#include "sc_utilities.h"

using namespace std;
using namespace SC;

/**
 * Customized test case for testing
 */
class HashTest : public ::testing::Test {
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
};

// Check how good hash_combine is
// You can test with many seeds at the same time
// The expected value of all seeds with the same m value
// is the same.
TEST_F(HashTest, test1) {
	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<float> real_dis(0.0, 100000.0);
	int seed1 = 0, seed2 = 0;
	float m = real_dis(gen);
	hash_combine<float>(seed1,m);
	hash_combine<float>(seed2,m);
	EXPECT_EQ(seed1,seed2);
	cout << seed1 << endl;
}

TEST_F(HashTest, test2) {
	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<float> real_dis(0.0, 100000.0);
	int seed1 = 0, seed2 = 0;
	float m1 = real_dis(gen);
	float m2 = real_dis(gen);
	EXPECT_NE(m1,m2);
	hash_combine<float>(seed1,m1);
	hash_combine<float>(seed2,m2);
	EXPECT_NE(seed1,seed2);
	cout << seed1 << endl;
	cout << seed2 << endl;
}

TEST_F(HashTest, test3) {
	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<ushort> dis(1,1 << 16 - 1);
	ushort u1[4], u2[4];
	for(int i = 0; i < 4; i++) {
		u1[i] = dis(gen);
		u2[i] = dis(gen);
		EXPECT_NE(u1[i],u2[i]);
	}
	int seed1 = vector_hash(&u1[0],4);
	int seed2 = vector_hash(&u2[0],4);
	EXPECT_NE(seed1,seed2);
	cout << seed1 << endl;
	cout << seed2 << endl;
}

TEST_F(HashTest, test4) {
	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<ushort> dis(1,1 << 16 - 1);
	ushort u[4];
	for(int i = 0; i < 4; i++) {
		u[i] = dis(gen);
	}
	int seed1 = vector_hash(&u[0],4);
	int seed2 = vector_hash(&u[0],4);
	EXPECT_EQ(seed1,seed2);
	cout << seed1 << endl;
}
