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

#pragma once

#include <string.h>
#include <iostream>

#include "filesystem/path.h"

namespace GALS
{
namespace UTILITIES
{
/*! \class FileUtils
 *
 * Class to handle file utilities, such as writing, reading.
 */
class FileUtils
{
 public:
  /*! Default constructor.
   */
  FileUtils();

  /*! Destructor
   */
  ~FileUtils();

  /*! Set root directory path.
   *
   * \param root_dir name of root directory.
   */
  void setRootDirectory(const std::string root_dir);

  /*! Get root directory path.
   *
   * \return root directory string.
   */
  const std::string getRootDirectory() const;

  /*! Checks if file exists or not.
   *
   * \param file_name name of file with relative or absolute path.
   *
   * \return true if exists, false otherwise.
   */
  bool fileExists(const std::string file_name) const;

  /*! Remove file.
   *
   * \param file_name name of file with relative or absolute path.
   *
   * \return true if success, false otherwise.
   */
  bool removeFile(const std::string file_name) const;

  /*! Create directory.
   *
   * This function creates any non-existing directories in the path.
   *
   * \param dir_name name of the directory to create.
   *
   * \return true if success, false otherwise.
   */
  bool createDirectory(const std::string dir_name) const;

  /*! Remove directory path.
   *
   * \param dir_name name of the path to remove.
   *
   * \return true if success, false otherwise.
   */
  bool removeDirectory(const std::string dir_name) const;

  /*! Write to file.
   *
   * \param file_name name of the file.
   * \param field variable whose data will be written to a file.
   */
  template <typename T>
  void write(const std::string file_name, const T& field);

 private:
  std::string m_root_dir;
};

}  // namespace UTILITIES
}  // namespace GALS
