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

#include <cpu/mat3.h>
#include <cpu/utilities.h>

#include <gtest/gtest.h>

template<typename T>
void test_subscript_operator(const GALS::CPU::Mat3<T> &mat3)
{
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         const T elem = mat3(i, j);
}

TEST(CPU, MAT3_INT)
{
   GALS::CPU::Mat3<int> mat3_1(1, 2, 3, 4, 5, 6, 7, 8, 9);
  GALS::CPU::Mat3<int> mat3;

  EXPECT_TRUE(mat3.size() == 9);

  mat3[0] = 9;
  EXPECT_TRUE(mat3[0] == 9);

  const int i = mat3[0];
  EXPECT_TRUE(i == 9);

  GALS::CPU::Mat3<int> mat_other;
  mat_other[0] = 1, mat_other[1] = 2, mat_other[2] = 3;

  // Assign mat_other to mat3.
  mat3 = mat_other;
  for (int i = 0; i < 9; ++i) EXPECT_TRUE(mat_other[i] == mat3[i]);

  mat_other(1, 1) = mat3(1, 1);

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) EXPECT_TRUE(mat_other(i, j) == mat3(i, j));

  EXPECT_TRUE(mat3 == mat_other);

  const int j = mat3(0, 2);
  EXPECT_TRUE(j == mat3(0, 2));

  test_subscript_operator<int>(mat3);
}

TEST(CPU, MAT3_DOUBLE)
{
   GALS::CPU::Mat3<double> mat3_1(1., 2., 3., 4, 5, 6, 7, 8, 9);

  GALS::CPU::Mat3<double> mat3;

  EXPECT_TRUE(mat3.size() == 9);

  mat3[0] = 9.;
  EXPECT_TRUE(GALS::is_equal(mat3[0], 9.));

  const double i = mat3[0];
  EXPECT_TRUE(GALS::is_equal(i, 9.));

  GALS::CPU::Mat3<double> mat_other;
  mat_other[0] = 1.1, mat_other[1] = 2.2, mat_other[2] = 3.3;

  // Assign mat_other to mat3.
  mat3 = mat_other;
  for (int i = 0; i < 9; ++i) EXPECT_TRUE(GALS::is_equal(mat_other[i], mat3[i]));

  mat_other(1, 1) = mat3(1, 1);

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) EXPECT_TRUE(GALS::is_equal(mat_other(i, j), mat3(i, j)));

  EXPECT_TRUE(mat3 == mat_other);

  const double j = mat3(0, 2);
  EXPECT_TRUE(GALS::is_equal(j, mat3(0, 2)));

  test_subscript_operator<double>(mat3);
}
