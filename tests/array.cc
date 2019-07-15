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

#include "gals/utilities/array.h"
#include "gals/utilities/utilities.h"
#include "gals/utilities/vec3.h"
#include "gals/utilities/vec_n.h"

#include <gtest/gtest.h>

#include <math.h>

// TODO (lakshman): All instances of array needs to have unit tests.
TEST(CPU, ARRAY)
{
  // scalar array on 1D grid.
  GALS::CPU::Grid<double, 1> grid(10, 1, 1);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, double> levelset(grid);

  EXPECT_TRUE(levelset.size() == 12);

  // Test grid().
  const GALS::CPU::Grid<double, 1> &grid_2 = levelset.grid();
  EXPECT_TRUE(grid_2.size() == 12);

  // 0D component array on 1D grid.
  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::VecN<double, 0>> twod_array_vecn0(grid);
  EXPECT_TRUE(twod_array_vecn0.size() == 12);

  // 2D component array on 1D grid.
  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::VecN<double, 2>> twod_array(grid);
  EXPECT_TRUE(twod_array.size() == 12);

  GALS::CPU::VecN<double, 2> vec;
  vec[0] = 0.9, vec[1] = 0.12;
  twod_array(0, 0, 0) = vec;
  EXPECT_TRUE(fabs(twod_array(0, 0, 0)[0] - vec[0]) < 1e-10);
  EXPECT_TRUE(fabs(twod_array(0, 0, 0)[1] - vec[1]) < 1e-10);

  // Test operator(node_id).
  GALS::CPU::Vec3<int> node_id(0, 0, 0);
  EXPECT_TRUE(fabs(twod_array(node_id)[0] - vec[0]) < 1e-10);

  EXPECT_TRUE(fabs(twod_array[1][0] - vec[0]) < 1e-10);

  twod_array[1][0] = 0.95;
  EXPECT_TRUE(fabs(twod_array[1][0] - 0.95) < 1e-10);

  EXPECT_TRUE(twod_array.numCells()[0] == 10);

  // const subscript operator.
  auto access = [&](const GALS::CPU::Array<GALS::CPU::Grid<double, 1>, double> &ls) {
    EXPECT_TRUE(ls.size() == 12);
    EXPECT_TRUE(GALS::is_equal(ls[2], -2.1));
  };

  levelset[2] = -2.1;
  access(levelset);

  // const paranthesis operator given node_id.
  auto paranthesis_access = [&](const GALS::CPU::Array<GALS::CPU::Grid<double, 1>, double> &ls) {
    const GALS::CPU::Vec3<int> node_id(1, 0, 0);
    EXPECT_TRUE(GALS::is_equal(ls(node_id), -2.2));
  };

  levelset(GALS::CPU::Vec3<int>(1, 0, 0)) = -2.2;
  paranthesis_access(levelset);

  const auto &data_levelset = levelset.data();
  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, double> levelset_2(grid);
  levelset_2 = levelset;
  paranthesis_access(levelset_2);
}
