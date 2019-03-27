/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include "grid.h"
#include <fstream>
#include <iostream>

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
GALS::CPU::VecN<T, 3>& GALS::CPU::Grid<T, DIM>::x(const int i, const int j, const int k)
{
  return m_grid[this->getIndex(i, j, k)];
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
const std::vector<int> GALS::CPU::Grid<T, DIM>::getNumCells() const
{
  return std::vector<int>{m_nx, m_ny, m_nz};
}

template <typename T, int DIM>
const int GALS::CPU::Grid<T, DIM>::getPadding() const
{
  return m_pad;
}

template <typename T, int DIM>
const std::size_t GALS::CPU::Grid<T, DIM>::getIndex(const int i, const int j, const int k)
{
  const std::size_t idx = ((k + m_pad) * (m_nx + 2 * m_pad) * (m_ny + 2 * m_pad) * m_mask[2]) +
                          ((j + m_pad) * (m_nx + 2 * m_pad) * m_mask[1]) + ((i + m_pad) * m_mask[0]);

  return idx;
}

template <typename T, int DIM>
const GALS::CPU::VecN<T, 3> GALS::CPU::Grid<T, DIM>::dX() const
{
  return m_dx;
}

template <typename T, int DIM>
GALS::CPU::VecN<T, 3>& GALS::CPU::Grid<T, DIM>::operator()(const int i, const int j, const int k)
{
  return m_grid[this->getIndex(i, j, k)];
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

  VecN<T, 3> elem;

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
void GALS::CPU::Grid<T, DIM>::writeToFile(std::string file_Name, std::string dir_Name, bool show_padding)
{
  std::string full_File_Name = dir_Name + "/" + file_Name;
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
        VecN<T, 3> coordinate = this->x(i, j, k);
        output_file << this->getIndex(i, j, k) << '\t' << coordinate[0] << '\t' << coordinate[1] << '\t'
                    << coordinate[2] << std::endl;
      }
    }
  }
}

template class GALS::CPU::Grid<double, 1>;
template class GALS::CPU::Grid<double, 2>;
template class GALS::CPU::Grid<double, 3>;
