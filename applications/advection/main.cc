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

#include <iostream>
#include <string>

#include "gals/analytical-fields/velocity.h"
#include "gals/input-parser.h"
#include "gals/utilities/file-utils.h"
#include "gals/utilities/grid.h"
#include "gals/utilities/vec3.h"

int main(int argc, char **argv)
{
  std::cout << "Inside applications/advection" << std::endl;

  std::string inputs_file;
  if (argc == 1) {
    std::cout << "./<executable> <inputs_file>" << std::endl;
    exit(0);
  } else {
    inputs_file = std::string(argv[1]);

    // TODO (lakshman): Check if the file exists.
  }

  const int dim = 2;
  using T = double;
  using TV = GALS::CPU::Vec3<T>;
  using T_GRID = GALS::CPU::Grid<T, dim>;

  GALS::INPUT_FIELDS::InputFields input_fields;
  GALS::INPUT_PARSER::InputParser input_parser;

  input_parser.parse(inputs_file, &input_fields);

  // Construct grid.
  const auto &grid_inputs = *(input_fields.m_grid);
  T_GRID grid(grid_inputs.nx, grid_inputs.ny, grid_inputs.nz);
  grid.generate(grid_inputs.x_min, grid_inputs.x_max, grid_inputs.y_min, grid_inputs.y_max, grid_inputs.z_min,
                grid_inputs.z_max);

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

  // Variable to store grid positions.
  GALS::CPU::Array<T_GRID, GALS::CPU::Vec3<T>> positions(grid);

  for (int i = i_min; i < i_max; ++i)
    for (int j = j_min; j < j_max; ++j)
      for (int k = k_min; k < k_max; ++k) {
        positions(i, j, k) = grid(i, j, k);
      }

  // Construct velocity.
  const auto &velocity_inputs = *(input_fields.m_velocity);
  GALS::CPU::Array<T_GRID, GALS::CPU::Vec3<T>> velocity_field(grid);
  GALS::ANALYTICAL_FIELDS::Velocity<T_GRID, T> velocity(grid, velocity_inputs);
  velocity.compute(positions, velocity_field);

  // Write velocity to a file.
  GALS::UTILITIES::FileUtils file_utils;
  file_utils.setRootDirectory("tmp/advection/");
  file_utils.createDirectory(file_utils.getRootDirectory());
  file_utils.write(std::string(file_utils.getRootDirectory() + "velocity"), velocity_field);

  return 0;
}
