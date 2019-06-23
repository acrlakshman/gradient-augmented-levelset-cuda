///////////////////////////////////////////////////////////////////////////////
// Copyright 2019 Lakshman Anumolu, Raunak Bardia.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
// may be used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#include "gals/utilities/vec_n.h"
#include "gals/utilities/utilities.h"

#include <gtest/gtest.h>

#include <iostream>

template <typename T, int dim>
void test_subscript_operator(const GALS::CPU::VecN<T, dim> &vec_n)
{
  for (int i = 0; i < vec_n.size(); ++i) {
    const T elem = vec_n[i];
  }
}

TEST(CPU, VEC_N_INT_1)
{
  GALS::CPU::VecN<int, 1> vec_n;

  EXPECT_TRUE(vec_n.size() == 1);

  vec_n[0] = 9;
  EXPECT_TRUE(vec_n[0] == 9);

  const int i = vec_n[0];
  EXPECT_TRUE(i == 9);

  test_subscript_operator<int, 1>(vec_n);

  // Overloaded output operator.
  std::cout << "VEC_N_INT_1 (<<): " << vec_n << std::endl;
}

TEST(CPU, VEC_N_INT_2)
{
  GALS::CPU::VecN<int, 2> vec_n;

  EXPECT_TRUE(vec_n.size() == 2);

  vec_n[0] = 9;
  EXPECT_TRUE(vec_n[0] == 9);

  const int i = vec_n[0];
  EXPECT_TRUE(i == 9);

  test_subscript_operator<int, 2>(vec_n);

  // Overloaded output operator.
  std::cout << "VEC_N_INT_2 (<<): " << vec_n << std::endl;
}

TEST(CPU, VEC_N_INT_3)
{
  GALS::CPU::VecN<int, 3> vec_n;

  EXPECT_TRUE(vec_n.size() == 3);

  vec_n[0] = 9;
  EXPECT_TRUE(vec_n[0] == 9);

  const int i = vec_n[0];
  EXPECT_TRUE(i == 9);

  test_subscript_operator<int, 3>(vec_n);

  // Overloaded output operator.
  std::cout << "VEC_N_INT_3 (<<): " << vec_n << std::endl;
}

TEST(CPU, VEC_N_DOUBLE_0)
{
  GALS::CPU::VecN<double, 0> vec_n;

  EXPECT_TRUE(vec_n.size() == 0);

  // Overloaded output operator.
  std::cout << "VEC_N_DOUBLE_0 (<<): " << vec_n << std::endl;
}

TEST(CPU, VEC_N_DOUBLE_1)
{
  GALS::CPU::VecN<double, 1> vec_n;

  EXPECT_TRUE(vec_n.size() == 1);

  vec_n[0] = 9.1;
  EXPECT_TRUE(GALS::is_equal(vec_n[0], 9.1));

  const double i = vec_n[0];
  EXPECT_TRUE(GALS::is_equal(i, 9.1));

  test_subscript_operator<double, 1>(vec_n);

  // Overloaded output operator.
  std::cout << "VEC_N_DOUBLE_1 (<<): " << vec_n << std::endl;
}

TEST(CPU, VEC_N_DOUBLE_2)
{
  GALS::CPU::VecN<double, 2> vec_n;

  EXPECT_TRUE(vec_n.size() == 2);

  vec_n[0] = 9.1;
  EXPECT_TRUE(GALS::is_equal(vec_n[0], 9.1));

  const double i = vec_n[0];
  EXPECT_TRUE(GALS::is_equal(i, 9.1));

  test_subscript_operator<double, 2>(vec_n);

  // Overloaded output operator.
  std::cout << "VEC_N_DOUBLE_2 (<<): " << vec_n << std::endl;
}

TEST(CPU, VEC_N_DOUBLE_3)
{
  GALS::CPU::VecN<double, 3> vec_n;

  EXPECT_TRUE(vec_n.size() == 3);

  vec_n[0] = 9.1;
  EXPECT_TRUE(GALS::is_equal(vec_n[0], 9.1));

  const double i = vec_n[0];
  EXPECT_TRUE(GALS::is_equal(i, 9.1));

  test_subscript_operator<double, 3>(vec_n);

  // Overloaded output operator.
  std::cout << "VEC_N_DOUBLE_3 (<<): " << vec_n << std::endl;
}

TEST(CPU, VEC_N_DOUBLE_4)
{
  GALS::CPU::VecN<double, 4> vec_n;

  EXPECT_TRUE(vec_n.size() == 4);

  vec_n[0] = 9.1;
  EXPECT_TRUE(GALS::is_equal(vec_n[0], 9.1));

  const double i = vec_n[0];
  EXPECT_TRUE(GALS::is_equal(i, 9.1));

  test_subscript_operator<double, 4>(vec_n);

  // Overloaded output operator.
  std::cout << "VEC_N_DOUBLE_4 (<<): " << vec_n << std::endl;
}
