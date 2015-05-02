/*
 * test_queue.cpp
 *
 *  Created on: 2014/12/29
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cstring>
#include <queue>
#include <random>
#include <gtest/gtest.h>
#include <boost/heap/priority_queue.hpp>
#include <boost/heap/binomial_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include "sc_utilities.h"

using namespace SC;

/**
 * Struct for store id and values
 */
typedef struct {
	int id;
	float value;
} p2;

typedef struct {
	int id1, id2, id3, id4;
	float value;
} p4;

class comp2
{
  bool reverse;
public:
  comp2(const bool& revparam=false)
    {reverse=revparam;}
  bool operator() (const p2& lhs, const p2&rhs) const
  {
    if (reverse) return (lhs.value>rhs.value);
    else return (lhs.value<rhs.value);
  }
};

class comp4
{
  bool reverse;
public:
  comp4(const bool& revparam=false)
    {reverse=revparam;}
  bool operator() (const p4& lhs, const p4&rhs) const
  {
    if (reverse) return (lhs.value>rhs.value);
    else return (lhs.value<rhs.value);
  }
};

/**
 * Customized test case for testing
 */
class QueueTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		N = 1 << 21;
		data = (float *)::operator new(N * sizeof(float));
		data2 = (float *)::operator new(N * sizeof(float));
		data_p = (p2 *)::operator new(N * sizeof(p2));
		id = (int *)::operator new(N * sizeof(int));
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<float> real_dis(0.0f, 100000.0f);
		for(int i = 0; i < N; i++) {
			data[i] = real_dis(gen);
			data_p[i].id = i;
			data_p[i].value = data[i];
		}
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
	static float * data, * data2;
	static int * id;
	static p2 * data_p;
	static int N;
};

float * QueueTest::data;
float * QueueTest::data2;
int * QueueTest::id;
p2 * QueueTest::data_p;
int QueueTest::N;
std::priority_queue<p2,vector<p2>,comp2> heap(comp2(false));
boost::heap::priority_queue<p2,boost::heap::compare<comp2>> heap2(comp2(false));
boost::heap::binomial_heap<p2,boost::heap::compare<comp2>> heap3(comp2(false));
boost::heap::fibonacci_heap<p2,boost::heap::compare<comp2>> heap4(comp2(false));

TEST_F(QueueTest, push1) {
	for(int i = 0; i < N; i++) {
		heap.push(data_p[i]);
	}
}

TEST_F(QueueTest, push2) {
	int sum = 0;
	for(int i = 0; i < N; i++) {
		pq_insert(data2,id,data[i],i,sum,N,false);
	}
}

TEST_F(QueueTest, push3) {
	for(int i = 0; i < N; i++) {
		heap2.push(data_p[i]);
	}
}

TEST_F(QueueTest, push4) {
	for(int i = 0; i < N; i++) {
		heap3.push(data_p[i]);
	}
}

TEST_F(QueueTest, push5) {
	for(int i = 0; i < N; i++) {
		heap4.push(data_p[i]);
	}
}

TEST_F(QueueTest, pop1) {
	float  d;
	for(int i = 0; i < N; i++) {
		d = heap.top().value;
		heap.pop();
	}
}

TEST_F(QueueTest, pop2) {
	int sum = N;
	for(int i = 0; i < N; i++) {
		pq_pop(data2,id,sum,false);
	}
}

TEST_F(QueueTest, pop3) {
	float  d;
	for(int i = 0; i < N; i++) {
		d = heap2.top().value;
		heap2.pop();
	}
}

TEST_F(QueueTest, pop4) {
	float  d;
	for(int i = 0; i < N; i++) {
		d = heap3.top().value;
		heap3.pop();
	}
}

TEST_F(QueueTest, pop5) {
	float  d;
	for(int i = 0; i < N; i++) {
		d = heap4.top().value;
		heap4.pop();
	}
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

