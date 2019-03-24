/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include <cpu/vec_n.h>

#include <gtest/gtest.h>

using namespace std;

TEST(CPU, VEC_N)
{
  GALS::CPU::VecN<int, 2> vec_n;

  EXPECT_TRUE(vec_n.size() == 2);

  vec_n[0] = 9;
  EXPECT_TRUE(vec_n[0] == 9);

  const int i = vec_n[0];
  EXPECT_TRUE(i == 9);
}
