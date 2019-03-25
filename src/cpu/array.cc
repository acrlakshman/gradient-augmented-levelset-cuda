/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include "array.h"
#include "vec_n.h"

template <typename T_GRID, typename T_ARRAY>
GALS::CPU::Array<T_GRID, T_ARRAY>::Array(Grid<typename T_GRID::value_type, T_GRID::dim>& grid)
    : m_grid(grid),
      m_nx(grid.getNumCells()[0]),
      m_ny(grid.getNumCells()[1]),
      m_nz(grid.getNumCells()[2]),
      m_pad(grid.getPadding())
{
  m_data.resize(m_grid.size());

  for (int i = 0; i < m_data.size(); ++i) {
    m_data[i] = T_ARRAY();
  }
}

template <typename T_GRID, typename T_ARRAY>
GALS::CPU::Array<T_GRID, T_ARRAY>::~Array()
{
  m_data.clear();
  m_data.shrink_to_fit();
}

template <typename T_GRID, typename T_ARRAY>
const std::size_t GALS::CPU::Array<T_GRID, T_ARRAY>::size() const
{
  return m_data.size();
}

template <typename T_GRID, typename T_ARRAY>
const T_ARRAY& GALS::CPU::Array<T_GRID, T_ARRAY>::operator[](const std::size_t idx) const
{
  return m_data[idx];
}

template <typename T_GRID, typename T_ARRAY>
T_ARRAY& GALS::CPU::Array<T_GRID, T_ARRAY>::operator[](const std::size_t idx)
{
  return m_data[idx];
}

template <typename T_GRID, typename T_ARRAY>
const T_ARRAY GALS::CPU::Array<T_GRID, T_ARRAY>::operator()(const int i, const int j, const int k) const
{
  const std::size_t idx = m_grid.getIndex(i, j, k);
  return m_data[idx];
}

template <typename T_GRID, typename T_ARRAY>
T_ARRAY& GALS::CPU::Array<T_GRID, T_ARRAY>::operator()(const int i, const int j, const int k)
{
  const std::size_t idx = m_grid.getIndex(i, j, k);
  return m_data[idx];
}

template class GALS::CPU::Array<GALS::CPU::Grid<double, 1>, double>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 2>, double>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 3>, double>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::VecN<int, 1>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::VecN<int, 2>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::VecN<int, 3>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 2>, GALS::CPU::VecN<int, 1>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 2>, GALS::CPU::VecN<int, 2>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 2>, GALS::CPU::VecN<int, 3>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 3>, GALS::CPU::VecN<int, 1>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 3>, GALS::CPU::VecN<int, 2>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 3>, GALS::CPU::VecN<int, 3>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::VecN<double, 1>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::VecN<double, 2>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::VecN<double, 3>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 2>, GALS::CPU::VecN<double, 1>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 2>, GALS::CPU::VecN<double, 2>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 2>, GALS::CPU::VecN<double, 3>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 3>, GALS::CPU::VecN<double, 1>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 3>, GALS::CPU::VecN<double, 2>>;
template class GALS::CPU::Array<GALS::CPU::Grid<double, 3>, GALS::CPU::VecN<double, 3>>;
