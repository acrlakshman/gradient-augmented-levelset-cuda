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

#include "gals/analytical-fields/velocity.h"

#include "gals/input-parser.h"
#include "gals/utilities/array.h"
#include "gals/utilities/file-utils.h"
#include "gals/utilities/utilities.h"
#include "gals/utilities/vec3.h"

#include <gtest/gtest.h>

#include <iostream>

namespace GU = GALS::UTILITIES;

/* * * * * *  TEST #1  * * * * * */
TEST(GALS, ANALYTICAL_FIELDS_VELOCITY_1D)
{
  typedef GALS::CPU::Grid<double, 1> T_GRID;

  // Initializing 1-D test grid.
  GALS::CPU::Grid<double, 1> grid(10, 1, 1);

  // grid generation
  grid.generate(-1, 1, -1, 1, -1, 1);

  // accessing grid details
  const auto mask = grid.getMask();
  const int pad = grid.getPadding();
  const auto num_cells = grid.numCells();

  // TODO incomplete.
}

/* * * * * *  TEST #2  * * * * * */
TEST(GALS, ANALYTICAL_FIELDS_VELOCITY_2D)
{
  using T = double;
  using T_GRID = GALS::CPU::Grid<T, 2>;

  // Initializing 2-D test grid.
  T_GRID grid(10, 10, 1);

  // grid generation.
  grid.generate(-1, 1, -1, 1, -1, 1);

  // accessing grid details
  const auto mask = grid.getMask();
  const int pad = grid.getPadding();
  const auto num_cells = grid.numCells();

  // defining working+ghost domain extent
  int i_min = -pad * mask[0];
  int j_min = -pad * mask[1];
  int k_min = -pad * mask[2];
  int i_max = num_cells[0] + pad * mask[0];
  int j_max = num_cells[1] + pad * mask[1];
  int k_max = num_cells[2] + pad * mask[2];

  // Input fields.
  GALS::INPUT_FIELDS::InputFields input_fields;

  GALS::INPUT_PARSER::InputParser input_parser;
  input_parser.parse("../../tests/inputs", &input_fields);

  const auto &velocity_inputs = *(input_fields.m_velocity);

  GALS::CPU::Array<T_GRID, GALS::CPU::Vec3<T>> positions(grid);

  for (int i = i_min; i < i_max; ++i)
    for (int j = j_min; j < j_max; ++j)
      for (int k = k_min; k < k_max; ++k) {
        positions(i, j, k) = grid(i, j, k);
      }

  GALS::CPU::Array<T_GRID, GALS::CPU::Vec3<T>> velocity_field(grid);
  GALS::ANALYTICAL_FIELDS::Velocity<T_GRID, T> velocity(grid, velocity_inputs);
  velocity.compute(positions, velocity_field);

  // TODO (lakshman): Compare these values with results from matlab.
  // for (int i = 0; i < num_cells[0]; ++i)
  // for (int j = 0; j < num_cells[1]; ++j)
  // for (int k = 0; k < num_cells[2]; ++k)
  // std::cout << "velocity(" << GALS::CPU::Vec3<int>(i, j, k) << "): " << velocity_field(i, j, k) << std::endl;

  // Write velocity to a file.
  GU::FileUtils file_utils;
  file_utils.setRootDirectory("tmp/velocity/");
  file_utils.createDirectory(file_utils.getRootDirectory());
  file_utils.write(std::string(file_utils.getRootDirectory() + "velocity"), velocity_field);
}
