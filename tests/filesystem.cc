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

#include "filesystem/path.h"

#include <gtest/gtest.h>

#include <math.h>
#include <fstream>

TEST(GALS, FILESYSTEM)
{
  namespace fs = filesystem;

  // Simple test.
  fs::path path_src("../../src");

  EXPECT_TRUE(path_src.is_directory());
  EXPECT_FALSE(path_src.is_file());

  // Create a temporary directory.
  fs::path path_tmp("./tmp");
  EXPECT_TRUE(fs::create_directory(path_tmp));

  // Open a file and write to it.
  std::string test_file_1 = "./tmp/test_file_1.txt";
  std::ofstream tmp_file_1(test_file_1);
  tmp_file_1 << "GALS cuda" << std::endl;
  tmp_file_1.close();

  path_src.set(test_file_1);

  // Remove file and directory.
  EXPECT_TRUE(path_src.remove_file());
  EXPECT_TRUE(fs::remove_directory(path_tmp));
}
