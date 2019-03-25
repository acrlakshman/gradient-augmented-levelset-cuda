/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include <cpu/array.h>
#include <cpu/vec_n.h>

#include <gtest/gtest.h>

#include <math.h>

TEST(CPU, ARRAY)
{
  // scalar array on 1D grid.
  GALS::CPU::Grid<double, 1> grid(10, 1, 1);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, double> levelset(grid);

  EXPECT_TRUE(levelset.size() == 12);

  // 2D component array on 1D grid.
  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::VecN<double, 2>> twod_array(grid);
  EXPECT_TRUE(twod_array.size() == 12);

  GALS::CPU::VecN<double, 2> vec;
  vec[0] = 0.9, vec[1] = 0.12;
  twod_array(0, 0, 0) = vec;
  EXPECT_TRUE(fabs(twod_array(0, 0, 0)[0] - vec[0]) < 1e-10);
  EXPECT_TRUE(fabs(twod_array(0, 0, 0)[1] - vec[1]) < 1e-10);

  EXPECT_TRUE(fabs(twod_array[1][0] - vec[0]) < 1e-10);

  twod_array[1][0] = 0.95;
  EXPECT_TRUE(fabs(twod_array[1][0] - 0.95) < 1e-10);
}
