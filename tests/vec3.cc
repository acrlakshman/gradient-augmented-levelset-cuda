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

#include "gals/utilities/vec3.h"
#include "gals/utilities/utilities.h"

#include <gtest/gtest.h>

#include <iostream>

TEST(CPU, VEC3_INT)
{
  GALS::CPU::Vec3<int> vec3;

  EXPECT_TRUE(vec3.size() == 3);

  vec3[0] = -1;
  vec3[1] = 1;
  vec3[2] = 2;
  EXPECT_TRUE(vec3.min() == -1);

  EXPECT_TRUE(GALS::is_equal(vec3.mag(), 2.449489742783178));

  vec3[0] = 9;
  EXPECT_TRUE(vec3[0] == 9);

  const int i = vec3[0];
  EXPECT_TRUE(i == 9);

  GALS::CPU::Vec3<int> vec_other;
  vec_other[0] = 1, vec_other[1] = 2, vec_other[2] = 3;

  // Assign vec_other to vec3.
  vec3 = vec_other;
  for (int i = 0; i < 3; ++i) EXPECT_TRUE(vec_other[i] == vec3[i]);

  EXPECT_TRUE(vec3 == vec_other);

  // Constructor with 3 input arguments.
  GALS::CPU::Vec3<int> vec_3(1, 1, 4);
  EXPECT_TRUE(GALS::is_equal(vec_3[0], 1));

  // Overloaded `-` operator.
  GALS::CPU::Vec3<int> vec3_1 = vec_3 - vec3;
  EXPECT_TRUE(GALS::is_equal(vec3_1[0], 0) && GALS::is_equal(vec3_1[1], -1) && GALS::is_equal(vec3_1[2], 1));

  // Overloaded `*` operator.
  GALS::CPU::Vec3<int> vec3_mult = vec_3 * vec3;
  EXPECT_TRUE(GALS::is_equal(vec3_mult[0], 1 * 1) && GALS::is_equal(vec3_mult[1], 1 * 2) &&
              GALS::is_equal(vec3_mult[2], 4 * 3));

  vec3_mult = vec_3 * 2;
  EXPECT_TRUE(GALS::is_equal(vec3_mult[0], 1 * 2) && GALS::is_equal(vec3_mult[1], 1 * 2) &&
              GALS::is_equal(vec3_mult[2], 4 * 2));

  // Overloaded `/` operator.
  GALS::CPU::Vec3<int> vec3_div = vec_3 / vec3;
  EXPECT_TRUE(GALS::is_equal(vec3_div[0], 1 / 1) && GALS::is_equal(vec3_div[1], 1 / 2) &&
              GALS::is_equal(vec3_div[2], 4 / 3));

  // Overloaded output operator.
  std::cout << "VEC3_INT (<<): " << vec3 << std::endl;

  // Instantiate vec3 object using standard vector.
  GALS::CPU::Vec3<int> vec3_using_vector(std::vector<int>{0, 1, 2});
  EXPECT_TRUE(GALS::is_equal(vec3_using_vector[0], 0) && GALS::is_equal(vec3_using_vector[1], 1) &&
              GALS::is_equal(vec3_using_vector[2], 2));
}

TEST(CPU, VEC3_DOUBLE)
{
  GALS::CPU::Vec3<double> vec3;

  EXPECT_TRUE(vec3.size() == 3);

  vec3[0] = -1.9;
  vec3[1] = 1.2;
  vec3[2] = 2.1;
  EXPECT_TRUE(GALS::is_equal(vec3.min(), -1.9));

  EXPECT_TRUE(GALS::is_equal(vec3.mag(), 3.0757112998459397));

  vec3[0] = 9.;
  EXPECT_TRUE(GALS::is_equal(vec3[0], 9.));

  const double i = vec3[0];
  EXPECT_TRUE(GALS::is_equal(i, 9.));

  GALS::CPU::Vec3<double> vec_other;
  vec_other[0] = 1.1, vec_other[1] = 2.2, vec_other[2] = 3.3;

  // Assign vec_other to vec3.
  vec3 = vec_other;
  for (int i = 0; i < 3; ++i) EXPECT_TRUE(GALS::is_equal(vec_other[i], vec3[i]));

  EXPECT_TRUE(vec3 == vec_other);

  // Constructor with 3 input arguments.
  GALS::CPU::Vec3<double> vec_3(1.1, 2.1, 3.1);
  EXPECT_TRUE(GALS::is_equal(vec_3[0], 1.1));

  // Overloaded `-` operator.
  GALS::CPU::Vec3<double> vec3_1 = vec_3 - vec3;
  EXPECT_TRUE(GALS::is_equal(vec3_1[0], 0.) && GALS::is_equal(vec3_1[1], -0.1) && GALS::is_equal(vec3_1[2], -0.2));

  // Overloaded `*` operator.
  GALS::CPU::Vec3<double> vec3_mult = vec_3 * vec3;
  EXPECT_TRUE(GALS::is_equal(vec3_mult[0], 1.1 * 1.1) && GALS::is_equal(vec3_mult[1], 2.1 * 2.2) &&
              GALS::is_equal(vec3_mult[2], 3.1 * 3.3));

  vec3_mult = vec_3 * 2.;
  EXPECT_TRUE(GALS::is_equal(vec3_mult[0], 1.1 * 2) && GALS::is_equal(vec3_mult[1], 2.1 * 2) &&
              GALS::is_equal(vec3_mult[2], 3.1 * 2));

  // Overloaded output operator.
  // Overloaded `/` operator.
  GALS::CPU::Vec3<double> vec3_div = vec_3 / vec3;
  EXPECT_TRUE(GALS::is_equal(vec3_div[0], 1.1 / 1.1) && GALS::is_equal(vec3_div[1], 2.1 / 2.2) &&
              GALS::is_equal(vec3_div[2], 3.1 / 3.3));

  // Overloaded output operator.
  std::cout << "VEC3_DOUBLE (<<): " << vec3 << std::endl;

  // Instantiate vec3 object using standard vector.
  GALS::CPU::Vec3<double> vec3_using_vector(std::vector<double>{0.2, 1.34, 2.12});
  EXPECT_TRUE(GALS::is_equal(vec3_using_vector[0], 0.2) && GALS::is_equal(vec3_using_vector[1], 1.34) &&
              GALS::is_equal(vec3_using_vector[2], 2.12));
}
