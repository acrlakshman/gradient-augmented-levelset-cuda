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

#include <fstream>

#include "gals/utilities/array.h"
#include "gals/utilities/grid.h"
#include "gals/utilities/utilities.h"
#include "gals/utilities/vec3.h"

namespace fs = filesystem;

GALS::UTILITIES::FileUtils::FileUtils() {}

GALS::UTILITIES::FileUtils::~FileUtils() {}

void GALS::UTILITIES::FileUtils::setRootDirectory(const std::string root_dir) { m_root_dir = root_dir; }

const std::string GALS::UTILITIES::FileUtils::getRootDirectory() const { return m_root_dir; }

bool GALS::UTILITIES::FileUtils::removeFile(const std::string file_name) const
{
  fs::path file_path(file_name);

  return file_path.remove_file();
}

bool GALS::UTILITIES::FileUtils::createDirectory(const std::string dir_name) const
{
  fs::path dir_path(dir_name);

  if (!dir_path.exists()) {
    return fs::create_directories(dir_path);
  }

  return true;
}

bool GALS::UTILITIES::FileUtils::removeDirectory(const std::string dir_name) const
{
  fs::path dir_path(dir_name);

  return fs::remove_directories(dir_path);
}

template <typename T>
void GALS::UTILITIES::FileUtils::write(const std::string file_name, const T& field)
{
  std::ofstream ofs(file_name);

  ofs << field;
  ofs.close();
}

#define P(...) __VA_ARGS__
#define _INSTANTIATE_WRITE_(type) template void GALS::UTILITIES::FileUtils::write<type>(const std::string, const type&);

_INSTANTIATE_WRITE_(P(GALS::CPU::Array<GALS::CPU::Grid<double, 1>, double>));
_INSTANTIATE_WRITE_(P(GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::Vec3<double>>));
