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
/*! \class Array
 *
 * Class to create Array.
 */
template <typename T_GRID, typename T_ARRAY>
class Array
{
 public:
  /*! Constructor called using grid.
   *
   * \param grid object of Grid.
   */
  Array(Grid<typename T_GRID::value_type, T_GRID::dim> &grid);

  /*! Destructor
   */
  ~Array();

  /*! Returns 1D array size of array.
   *
   * Data of array is stored in 1D vector with some padding. This function returns 1D vector size.
   *
   * \return 1D array size.
   */
  const std::size_t size() const;

  /*! Overloaded subscript operator to return value of array at a given 1D array based index.
   *
   * \param idx 1D array based index.
   *
   * \return value at given index.
   */
  const T_ARRAY &operator[](const std::size_t idx) const;

  /*! Overloaded subscript operator to return reference to value of array at a given 1D array based index.
   *
   * \param idx 1D array based index.
   *
   * \return reference to value at given index.
   */
  T_ARRAY &operator[](const std::size_t idx);

  /*! Overloaded operator to return value of array using 3D cell index.
   *
   * \param i zero based index along x-direction.
   * \param j zero based index along y-direction.
   * \param k zero based index along z-direction.
   *
   * \return value at given 3D cell index.
   */
  const T_ARRAY operator()(const int i, const int j, const int k) const;

  /*! Overloaded operator to return reference to value of array using 3D cell index.
   *
   * \param i zero based index along x-direction.
   * \param j zero based index along y-direction.
   * \param k zero based index along z-direction.
   *
   * \return reference to value at given 3D cell index.
   */
  T_ARRAY &operator()(const int i, const int j, const int k);

 private:
  T_GRID &m_grid;
  const int m_nx, m_ny, m_nz, m_pad;
  std::vector<T_ARRAY> m_data;
};

}  // namespace CPU
}  // namespace GALS
