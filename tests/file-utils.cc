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

#include "gals/utilities/file-utils.h"
#include "gals/utilities/array.h"
#include "gals/utilities/grid.h"
#include "gals/utilities/vec3.h"

#include <gtest/gtest.h>

namespace GU = GALS::UTILITIES;

TEST(GALS, UTILITIES_FILE_UTILS)
{
  GU::FileUtils file_utils;
  file_utils.setRootDirectory("root_dir/proj/case/");

  // Create project directory.
  EXPECT_TRUE(file_utils.createDirectory(file_utils.getRootDirectory()));

  // scalar array on 1D grid.
  GALS::CPU::Grid<double, 1> grid(10, 1, 1);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::Vec3<double>> oned_array_vec3(grid);

  // Write array to file.
  file_utils.write(std::string(file_utils.getRootDirectory() + "file.txt"), oned_array_vec3);

  // Remove file.
  EXPECT_TRUE(file_utils.removeFile(std::string(file_utils.getRootDirectory() + "file.txt")));

  // Remove entire directory.
  EXPECT_TRUE(file_utils.removeDirectory(file_utils.getRootDirectory()));
}
