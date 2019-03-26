/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#pragma once

#include "vec_n.h"

#include <vector>
#include <string>

namespace GALS
{
namespace CPU
{
/*! \class Grid
 *
 * Class to create grid.
 */
template <typename T, int DIM = 3>
class Grid
{
 public:
  typedef T value_type;
  static const int dim = DIM;

  /*! Constructor with arguments for number of cells across x, y, z.
   *
   * \param nx No. of cells across x-direction.
   * \param ny No. of cells across y-direction.
   * \param nz No. of cells across z-direction.
   */
  Grid(int nx, int ny, int nz);

  /*! Constructor with arguments for number of cells across x, y.
   *
   * Number of cells across z-direction is defaulted to 1.
   *
   * \param nx No. of cells across x-direction.
   * \param ny No. of cells across y-direction.
   */
  Grid(int nx, int ny);

  /*! Constructor with arguments for number of cells across x.
   *
   * Number of cells across y, z-directions are defaulted to 1.
   *
   * \param nx No. of cells across x-direction.
   */
  Grid(int nx);

  /*! Destructor
   */
  ~Grid();

  /*! Returns position at a given index.
   *
   * \param i zero based index along x-direction.
   * \param j zero based index along y-direction.
   * \param k zero based index along z-direction.
   *
   * \return 3D position vector.
   */
  VecN<T, 3>& x(const int i, const int j = 0, const int k = 0);

  /*! Returns dimension of grid.
   *
   * \return dimension of grid.
   */
  const int dimension() const;

  /*! Returns 1D array size of grid.
   *
   * Data of grid is stored in 1D array with some padding. This function returns 1D array size.
   *
   * \return 1D array size.
   */
  const int size() const;

  /*! Returns pointer to mask.
   *
   * For efficient computation of 3D index to 1D index a mask is used to differentiate between 1, 2, 3 dimentions.
   * - mask = {1, 0, 0} for 1D
   * - mask = {1, 1, 0} for 2D
   * - mask = {1, 1, 1} for 3D
   *
   * \return integer pointer that points to integer array of size 3.
   */
  const int* getMask() const;

  /*! Returns a vector of size 3 with number of cells along x, y, z directions.
   *
   * Vector of size 3: {nx, ny, nz}.
   *
   * \return vector of size 3.
   */
  const std::vector<int> getNumCells() const;

  /*! Return current padding.
   *
   * \return padding.
   */
  const int getPadding() const;

  /*! Returns 1D index in stored array.
   *
   * 3D to 1D index mapping.
   *
   * \return 1D index.
   */
  const std::size_t getIndex(const int i, const int j, const int k);

  /*! Returns cell size which is a vector of size 3.
   *
   * \return cell size.
   */
  const VecN<T, 3> dX() const;

  /*! Operator overloaded to return co-ordinate values at a given 3D index.
   *
   * \return position.
   */
  VecN<T, 3>& operator()(const int i, const int j, const int k);

  /*! Set new padding value.
   *
   * When used this, `generate(...)` function must be called immediately after setting new padding.
   *
   * \param pad new padding value.
   */
  void setPadding(const int pad);

  /*! Generate grid using domain bounding box.
   *
   * \param x_min minimum x-coordinate.
   * \param x_max maximum x-coordinate.
   * \param y_min minimum y-coordinate.
   * \param y_max maximum y-coordinate.
   * \param z_min minimum z-coordinate.
   * \param z_max maximum z-coordinate.
   */
  void generate(T x_min, T x_max, T y_min, T y_max, T z_min, T z_max);

  /*! Print grid for debugging.
   *
   * \param file_name writes grid file with the prescribed name.
   * \param dir_name writes grid in the prescribed directory.
   * \param show_padding writes grid with or without padding cells.
   */
  void writeToFile(std::string file_name = "grid.dat", std::string dir_name = ".", bool show_padding = false);

 private:
  int m_dimension, m_nx, m_ny, m_nz, m_pad;
  int m_mask[3];
  std::vector<VecN<T, 3>> m_grid;
  VecN<T, 3> m_dx;
};

}  // namespace CPU
}  // namespace GALS
