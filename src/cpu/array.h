/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#pragma once

#include "grid.h"

#include <vector>

namespace GALS
{
namespace CPU
{
template <typename T_GRID, typename T_ARRAY>
class Array
{
 public:
  Array(Grid<typename T_GRID::value_type, T_GRID::dim> &grid);

  ~Array();

  const std::size_t size() const;

  const T_ARRAY &operator[](const std::size_t idx) const;

  T_ARRAY &operator[](const std::size_t idx);

  const T_ARRAY operator()(const int i, const int j, const int k) const;

  T_ARRAY &operator()(const int i, const int j, const int k);

 private:
  T_GRID &m_grid;
  const int m_nx, m_ny, m_nz, m_pad;
  std::vector<T_ARRAY> m_data;
};

}  // namespace CPU
}  // namespace GALS
