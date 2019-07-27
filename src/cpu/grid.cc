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

#include "gals/utilities/grid.h"

#include <limits.h>
#include <math.h>
#include <fstream>

#include "gals/input-parser.h"

template <typename T, int DIM>
GALS::CPU::Grid<T, DIM>::Grid(int nx, int ny, int nz)
    : m_dimension(DIM), m_nx(nx), m_ny(ny), m_nz(nz), m_pad(1), m_total_cells(nx * ny * nz)
{
  m_mask[0] = 1, m_mask[1] = 0, m_mask[2] = 0;

  if (dim == 2)
    m_mask[1] = 1;
  else if (dim == 3)
    m_mask[1] = 1, m_mask[2] = 1;

  m_box_min = Vec3<T>(INT_MIN, INT_MIN, INT_MIN);
  m_box_max = Vec3<T>(INT_MAX, INT_MAX, INT_MAX);

  m_dx = GALS::CPU::Vec3<T>(INT_MAX, INT_MAX, INT_MAX), m_one_over_dx = GALS::CPU::Vec3<T>(INT_MIN, INT_MIN, INT_MIN);
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
const GALS::CPU::Vec3<int> GALS::CPU::Grid<T, DIM>::numCells() const
{
  return Vec3<int>(m_nx, m_ny, m_nz);
}

template <typename T, int DIM>
const size_t GALS::CPU::Grid<T, DIM>::totalCells() const
{
  return m_total_cells;
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
const std::size_t GALS::CPU::Grid<T, DIM>::index(const Vec3<int> node_id) const
{
  return this->index(node_id[0], node_id[1], node_id[2]);
}

template <typename T, int DIM>
const GALS::CPU::Vec3<int> GALS::CPU::Grid<T, DIM>::baseNodeId(const Vec3<T>& x) const
{
  Vec3<int> base_node_id(INT_MAX, INT_MAX, INT_MAX);

  for (int axis = 0; axis < DIM; ++axis) {
    base_node_id[axis] = floor(((x[axis] - m_box_min[axis]) * m_one_over_dx[axis]) - 0.5);
  }

  return base_node_id;
}

template <typename T, int DIM>
const GALS::CPU::Vec3<T>& GALS::CPU::Grid<T, DIM>::dX() const
{
  return m_dx;
}

template <typename T, int DIM>
const GALS::CPU::Vec3<T>& GALS::CPU::Grid<T, DIM>::oneOverDX() const
{
  return m_one_over_dx;
}

template <typename T, int DIM>
const GALS::CPU::Vec3<T>& GALS::CPU::Grid<T, DIM>::operator()(const int i, const int j, const int k) const
{
  return m_grid[this->index(i, j, k)];
}

template <typename T, int DIM>
const GALS::CPU::Vec3<T>& GALS::CPU::Grid<T, DIM>::operator()(const Vec3<int> node_id) const
{
  return m_grid[this->index(node_id)];
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

  m_box_min[0] = x_min, m_box_min[1] = y_min, m_box_min[2] = z_min;
  m_box_max[0] = x_max, m_box_max[1] = y_max, m_box_max[2] = z_max;

  m_grid.resize((m_nz + 2 * m_pad * m_mask[2]) * (m_ny + 2 * m_pad * m_mask[1]) * (m_nx + 2 * m_pad * m_mask[0]));

  std::vector<T> domain_min({x_min, y_min, z_min}), domain_min_new(3);

  m_dx[0] = (x_max - x_min) / m_nx;
  m_dx[1] = (y_max - y_min) / m_ny;
  m_dx[2] = (z_max - z_min) / m_nz;

  for (int i = 0; i < DIM; ++i) m_one_over_dx[i] = static_cast<T>(1.) / m_dx[i];

  for (int i = 0; i < DIM; ++i) domain_min_new[i] = domain_min[i] + (m_dx[i] * 0.5) - (m_dx[i] * m_mask[i]);

  Vec3<T> elem;

  int i_min = -m_pad * m_mask[0], j_min = -m_pad * m_mask[1], k_min = -m_pad * m_mask[2];

  GALS::CPU::Vec3<int> vec;

  for (int i = i_min; i < m_nx + m_pad * m_mask[0]; ++i) {
    for (int j = j_min; j < m_ny + m_pad * m_mask[1]; ++j) {
      for (int k = k_min; k < m_nz + m_pad * m_mask[2]; ++k) {
        vec[0] = i - i_min;
        vec[1] = j - j_min;
        vec[2] = k - k_min;

        for (int axis = 0; axis < DIM; ++axis) elem[axis] = domain_min_new[axis] + vec[axis] * m_dx[axis];

        this->x(i, j, k) = elem;
      }
    }
  }

  // Update total cells.
  m_total_cells = 1;
  for (int d = 0; d < DIM; ++d) m_total_cells *= this->numCells()[d];
}

template class GALS::CPU::Grid<double, 1>;
template class GALS::CPU::Grid<double, 2>;
template class GALS::CPU::Grid<double, 3>;
