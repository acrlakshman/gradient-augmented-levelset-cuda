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

#include "grid.h"
#include "input-parser.h"

#include <fstream>

template <typename T, int DIM>
GALS::CPU::Grid<T, DIM>::Grid(int nx, int ny, int nz) : m_dimension(DIM), m_nx(nx), m_ny(ny), m_nz(nz), m_pad(1)
{
  m_mask[0] = 1, m_mask[1] = 0, m_mask[2] = 0;

  if (dim == 2)
    m_mask[1] = 1;
  else if (dim == 3)
    m_mask[1] = 1, m_mask[2] = 1;
}

template <typename T, int DIM>
GALS::CPU::Grid<T, DIM>::Grid(int nx, int ny) : Grid(nx, ny, 1)
{
}

template <typename T, int DIM>
GALS::CPU::Grid<T, DIM>::Grid(int nx) : Grid(nx, 1, 1)
{
}

template <typename T, int DIM>
GALS::CPU::Grid<T, DIM>::~Grid()
{
  m_grid.clear();
  m_grid.shrink_to_fit();
}

template <typename T, int DIM>
GALS::CPU::Vec3<T>& GALS::CPU::Grid<T, DIM>::x(const int i, const int j, const int k)
{
  return m_grid[this->index(i, j, k)];
}

template <typename T, int DIM>
const int GALS::CPU::Grid<T, DIM>::dimension() const
{
  return m_dimension;
}

template <typename T, int DIM>
const int GALS::CPU::Grid<T, DIM>::size() const
{
  return m_grid.size();
}

template <typename T, int DIM>
const int* GALS::CPU::Grid<T, DIM>::getMask() const
{
  return m_mask;
}

template <typename T, int DIM>
const std::vector<int> GALS::CPU::Grid<T, DIM>::numCells() const
{
  return std::vector<int>{m_nx, m_ny, m_nz};
}

template <typename T, int DIM>
const int GALS::CPU::Grid<T, DIM>::getPadding() const
{
  return m_pad;
}

template <typename T, int DIM>
const std::size_t GALS::CPU::Grid<T, DIM>::index(const int i, const int j, const int k) const
{
  const std::size_t idx = ((k + m_pad) * (m_nx + 2 * m_pad) * (m_ny + 2 * m_pad) * m_mask[2]) +
                          ((j + m_pad) * (m_nx + 2 * m_pad) * m_mask[1]) + ((i + m_pad) * m_mask[0]);

  return idx;
}

template <typename T, int DIM>
const GALS::CPU::Vec3<T>& GALS::CPU::Grid<T, DIM>::dX() const
{
  return m_dx;
}

template <typename T, int DIM>
const GALS::CPU::Vec3<T>& GALS::CPU::Grid<T, DIM>::operator()(const int i, const int j, const int k) const
{
  return m_grid[this->index(i, j, k)];
}

template <typename T, int DIM>
void GALS::CPU::Grid<T, DIM>::setPadding(const int pad)
{
  m_pad = pad;
}

template <typename T, int DIM>
void GALS::CPU::Grid<T, DIM>::generate(T x_min, T x_max, T y_min, T y_max, T z_min, T z_max)
{
  if (m_grid.size()) m_grid.clear(), m_grid.shrink_to_fit();

  m_grid.resize((m_nz + 2 * m_pad * m_mask[2]) * (m_ny + 2 * m_pad * m_mask[1]) * (m_nx + 2 * m_pad * m_mask[0]));

  std::vector<T> domain_min({x_min, y_min, z_min}), domain_min_new(3);

  m_dx[0] = (x_max - x_min) / m_nx;
  m_dx[1] = (y_max - y_min) / m_ny;
  m_dx[2] = (z_max - z_min) / m_nz;

  for (int i = 0; i < 3; ++i) domain_min_new[i] = domain_min[i] + (m_dx[i] * 0.5) - (m_dx[i] * m_mask[i]);

  Vec3<T> elem;

  int i_min = -m_pad * m_mask[0], j_min = -m_pad * m_mask[1], k_min = -m_pad * m_mask[2];

  for (int i = i_min; i < m_nx + m_pad * m_mask[0]; ++i) {
    for (int j = j_min; j < m_ny + m_pad * m_mask[1]; ++j) {
      for (int k = k_min; k < m_nz + m_pad * m_mask[2]; ++k) {
        elem[0] = domain_min_new[0] + (i - i_min) * m_dx[0];
        elem[1] = domain_min_new[1] + (j - j_min) * m_dx[1];
        elem[2] = domain_min_new[2] + (k - k_min) * m_dx[2];

        this->x(i, j, k) = elem;
      }
    }
  }
}

template <typename T, int DIM>
void GALS::CPU::Grid<T, DIM>::writeToFile(std::string file_name, std::string dir_name, bool show_padding)
{
  std::string full_File_Name = dir_name + "/" + file_name;
  std::ofstream output_file(full_File_Name);
  output_file << "Node_ID\tX\tY\tZ\n";

  int i_min = show_padding ? -m_pad * m_mask[0] : 0;
  int j_min = show_padding ? -m_pad * m_mask[1] : 0;
  int k_min = show_padding ? -m_pad * m_mask[2] : 0;
  int i_max = show_padding ? m_nx + m_pad * m_mask[0] : m_nx;
  int j_max = show_padding ? m_ny + m_pad * m_mask[1] : m_ny;
  int k_max = show_padding ? m_nz + m_pad * m_mask[2] : m_nz;

  for (int i = i_min; i < i_max; ++i) {
    for (int j = j_min; j < j_max; ++j) {
      for (int k = k_min; k < k_max; ++k) {
        Vec3<T> coordinate = this->x(i, j, k);
        output_file << this->index(i, j, k) << '\t' << coordinate[0] << '\t' << coordinate[1] << '\t' << coordinate[2]
                    << std::endl;
      }
    }
  }
}

template class GALS::CPU::Grid<double, 1>;
template class GALS::CPU::Grid<double, 2>;
template class GALS::CPU::Grid<double, 3>;
