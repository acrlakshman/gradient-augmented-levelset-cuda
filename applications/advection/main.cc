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

#include <algorithm>
#include <iostream>
#include <string>

#include "gals/analytical-fields/levelset.h"
#include "gals/analytical-fields/velocity.h"
#include "gals/cpu/gradient.h"
#include "gals/cpu/levelset.h"
#include "gals/cpu/temporal.h"
#include "gals/input-parser.h"
#include "gals/utilities/file-utils.h"
#include "gals/utilities/grid.h"
#include "gals/utilities/vec3.h"

const int dim = 2;
using T = double;
using TV = GALS::CPU::Vec3<T>;
using T_GRID = GALS::CPU::Grid<T, dim>;

static T compute_dt(const T cfl_max, const T_GRID &grid, const GALS::CPU::Array<T_GRID, TV> &velocity);

int main(int argc, char **argv)
{
  std::cout << "Inside applications/advection" << std::endl;

  std::string inputs_file;
  if (argc == 1) {
    std::cout << "<path-to-executable>/<executable> <path-to-inputs_file>" << std::endl;
    exit(0);
  } else {
    inputs_file = std::string(argv[1]);

    GALS::UTILITIES::FileUtils file_utils;
    if (!file_utils.fileExists(inputs_file)) {
      std::cout << "File: " << inputs_file << " does not exist" << std::endl;
      exit(0);
    }
  }

  GALS::INPUT_FIELDS::InputFields input_fields;
  GALS::INPUT_PARSER::InputParser input_parser;

  input_parser.parse(inputs_file, &input_fields);

  GALS::UTILITIES::FileUtils file_utils;

  const auto &general_inputs = *(input_fields.m_general);
  file_utils.setRootDirectory(general_inputs.output_directory + "/");

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
  GALS::CPU::LevelsetVelocity<T_GRID, T> levelset_velocity(grid);
  const auto &velocity_field = levelset_velocity.velocity();
  GALS::ANALYTICAL_FIELDS::Velocity<T_GRID, T> analytical_velocity(grid, velocity_inputs);

  // Build levelset.
  const auto &levelset_inputs = *(input_fields.m_levelset);
  GALS::CPU::Levelset<T_GRID, T> levelset(grid);
  auto &phi = levelset.phi();
  auto &psi = levelset.psi();
  auto &phi_prev = levelset.phiPrev();
  auto &psi_prev = levelset.psiPrev();
  auto &phi_mixed_derivatives_prev = levelset.phiMixedDerivativesPrev();
  GALS::ANALYTICAL_FIELDS::Levelset<T_GRID, T> analytical_levelset(grid, levelset_inputs);
  auto iii = 0;
  std::cout << "here " << ++iii << std::endl;

  analytical_levelset.compute(positions, levelset);
  std::cout << "here " << ++iii << std::endl;

  // Build time data.
  const auto &time_inputs = *(input_fields.m_time);

  T t_start = time_inputs.start;
  T t_end = time_inputs.end;
  T dt = time_inputs.dt;
  bool is_dt_fixed = std::strcmp(time_inputs.constant_dt.c_str(), "NO");
  T cfl_max = time_inputs.cfl_max;
  T next_write_time = time_inputs.write_interval;
  T sim_time = static_cast<T>(0);

  if (!is_dt_fixed) {
    dt = grid.dX().min() * static_cast<T>(0.5);
  }
  std::cout << "here " << ++iii << std::endl;

  // Compute gradient of levelset field.
  GALS::CPU::Gradient<T, T_GRID, GALS::CPU::ThirdOrder<T, T_GRID>>::compute(phi, psi);

  // Time loop
  int write_count = 1;
  bool run_sim = true;
  while (run_sim) {
    // Compute velocity and its gradient at current time.
    analytical_velocity.compute(positions, sim_time, levelset_velocity);

    if (!is_dt_fixed) dt = compute_dt(cfl_max, grid, velocity_field);
    std::cout << "=======================================\n" << std::endl;
    std::cout << "Time step = " << dt << std::endl;
    sim_time += dt;

    iii = 1000;
    std::cout << "here " << ++iii << std::endl;

    // Compute gradient of levelset field.
    GALS::CPU::Gradient<T, T_GRID, GALS::CPU::ThirdOrder<T, T_GRID>>::compute(phi, psi);
    phi_prev = phi;
    psi_prev = psi;
    // std::cout << "phi_prev = " << phi_prev << std::endl;
    // std::cout << "psi_prev = " << psi_prev << std::endl;
    std::cout << "here " << ++iii << std::endl;

    // Compute levelset mixed derivatives for `_prev` fields.
    levelset.computeMixedDerivatives(psi_prev, phi_mixed_derivatives_prev);
    // std::cout << "phi_mixed_derivatives_prev = " << phi_mixed_derivatives_prev << std::endl;
    std::cout << "here " << ++iii << std::endl;

    // Advect levelset.
    GALS::CPU::Temporal<T, T_GRID, GALS::TEMPORAL_SCHEMES::SEMI_LAGRANGIAN::Euler<T, T_GRID>>::compute(
        dt, levelset_velocity, levelset);
    std::cout << "here " << ++iii << std::endl;

    if (sim_time >= next_write_time) {
      auto write_dir = file_utils.getRootDirectory() + std::to_string(write_count);
      file_utils.createDirectory(write_dir);

      // Write velocity to a file.
      file_utils.write(std::string(write_dir + "/velocity"), velocity_field);

      // Write levelset to a file.
      file_utils.write(std::string(write_dir + "/phi"), levelset.phi());

      ++write_count;
      next_write_time += time_inputs.write_interval;
    }

    if (GALS::is_equal(sim_time, t_end) || sim_time > t_end) run_sim = false;
  }

  return 0;
}

T compute_dt(const T cfl_max, const T_GRID &grid, const GALS::CPU::Array<T_GRID, TV> &velocity)
{
  T dt = 0.;

  const auto dx = grid.dX();
  const auto dx_min = dx.min();

  // accessing grid details
  const auto mask = grid.getMask();
  const int pad = grid.getPadding();
  const auto num_cells = grid.numCells();

  // defining working+ghost domain extent
  int i_min = 0;
  int j_min = 0;
  int k_min = 0;
  int i_max = num_cells[0];
  int j_max = num_cells[1];
  int k_max = num_cells[2];

  T max_u_mag = static_cast<T>(0.);

  for (int i = i_min; i < i_max; ++i)
    for (int j = j_min; j < j_max; ++j)
      for (int k = k_min; k < k_max; ++k) {
        // std::cout << "velocity(" << i << "," << j << "," << k << "): " << velocity(i, j, k) << std::endl;
        auto vel = velocity(i, j, k).mag();
        max_u_mag = std::max(max_u_mag, vel);
      }
  T one_by_max_u_mag = static_cast<T>(1.) / max_u_mag;

  std::cout << "one_by_max_u_mag = " << one_by_max_u_mag << std::endl;
  std::cout << "\tcfl_max = " << cfl_max << std::endl;
  std::cout << "\tdx = " << dx << std::endl;
  std::cout << "\tdx_min = " << dx_min << std::endl;
  dt = cfl_max * dx_min * one_by_max_u_mag;

  return dt;
}
