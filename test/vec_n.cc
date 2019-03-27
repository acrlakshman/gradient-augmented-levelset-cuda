/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include <cpu/utilities.h>
#include <cpu/vec_n.h>

#include <gtest/gtest.h>

using namespace std;

TEST(CPU, VEC_N_INT_1)
{
  GALS::CPU::VecN<int, 1> vec_n;

  EXPECT_TRUE(vec_n.size() == 1);

  vec_n[0] = 9;
  EXPECT_TRUE(vec_n[0] == 9);

  const int i = vec_n[0];
  EXPECT_TRUE(i == 9);
}

TEST(CPU, VEC_N_INT_2)
{
  GALS::CPU::VecN<int, 2> vec_n;

  EXPECT_TRUE(vec_n.size() == 2);

  vec_n[0] = 9;
  EXPECT_TRUE(vec_n[0] == 9);

  const int i = vec_n[0];
  EXPECT_TRUE(i == 9);
}

TEST(CPU, VEC_N_INT_3)
{
  GALS::CPU::VecN<int, 3> vec_n;

  EXPECT_TRUE(vec_n.size() == 3);

  vec_n[0] = 9;
  EXPECT_TRUE(vec_n[0] == 9);

  const int i = vec_n[0];
  EXPECT_TRUE(i == 9);
}

TEST(CPU, VEC_N_DOUBLE_1)
{
  GALS::CPU::VecN<double, 1> vec_n;

  EXPECT_TRUE(vec_n.size() == 1);

  vec_n[0] = 9.1;
  EXPECT_TRUE(GALS::is_equal<double>(vec_n[0], 9.1));

  const double i = vec_n[0];
  EXPECT_TRUE(GALS::is_equal<double>(i, 9.1));
}

TEST(CPU, VEC_N_DOUBLE_2)
{
  GALS::CPU::VecN<double, 2> vec_n;

  EXPECT_TRUE(vec_n.size() == 2);

  vec_n[0] = 9.1;
  EXPECT_TRUE(GALS::is_equal<double>(vec_n[0], 9.1));

  const double i = vec_n[0];
  EXPECT_TRUE(GALS::is_equal<double>(i, 9.1));
}

TEST(CPU, VEC_N_DOUBLE_3)
{
  GALS::CPU::VecN<double, 3> vec_n;

  EXPECT_TRUE(vec_n.size() == 3);

  vec_n[0] = 9.1;
  EXPECT_TRUE(GALS::is_equal<double>(vec_n[0], 9.1));

  const double i = vec_n[0];
  EXPECT_TRUE(GALS::is_equal<double>(i, 9.1));
}
