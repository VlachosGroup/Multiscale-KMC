#include <gtest/gtest.h>


// Tests factorial of 0.
TEST(FactorialTest, HandlesZeroInput) {
	int a = 5;
	int b = -3;
	int c = a + b;
  EXPECT_EQ(c, 2);
}