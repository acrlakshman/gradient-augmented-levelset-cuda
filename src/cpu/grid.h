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

#include <string>
#include <vector>

#include "input-fields/grid.h"
#include "mat3.h"
#include "vec3.h"

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
  typedef Vec3<T> position_type;
  static const int dim = DIM;

  //! vector along x: {1, 0, 0}; y: {0, 1, 0}; z: {0, 0, 1}.
  static const Mat3<int> axis_vectors;

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
  Vec3<T>& x(const int i, const int j = 0, const int k = 0);

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
  const std::vector<int> numCells() const;

  /*! Returns total number of cells excluding ghost (padded) cells.
   *
   * \return number of cells of type size_t.
   */
  const size_t totalCells() const;

  /*! Return current padding.
   *
   * \return padding.
   */
  const int getPadding() const;

  /*! Returns 1D index in stored array.
   *
   * 3D to 1D index mapping.
   *
   * \param i zero based index along x-direction.
   * \param j zero based index along y-direction.
   * \param k zero based index along z-direction.
   *
   * \return 1D index.
   */
  const std::size_t index(const int i, const int j, const int k) const;

  /*! Given a node id, returns 1D index in stored array.
   *
   * \param node_id node id of type GALS::CPU::Vec3<int>.
   *
   * \return 1D index.
   */
  const std::size_t index(const Vec3<int> node_id) const;

  /*! Return base node id (i, j, k) that encloses given position.
   *
   * \param x position.
   *
   * \return base node id (i, j, k).
   */
  const Vec3<int> baseNodeId(const Vec3<T>& x) const;

  /*! Returns cell size which is a vector of size 3.
   *
   * \return cell size.
   */
  const Vec3<T>& dX() const;

  /*! Returns one over cell size which is a vector of size 3.
   *
   * \return one over cell size.
   */
  const Vec3<T>& oneOverDX() const;

  /*! Operator overloaded to return co-ordinate values at a given 3D index.
   *
   * \return position.
   */
  const Vec3<T>& operator()(const int i, const int j, const int k) const;

  /*! Operator overloaded to return co-ordinate values at a given 3D index.
   *
   * \param node_id node index of type NodeId.
   *
   * \return position.
   */
  const Vec3<T>& operator()(const Vec3<int> node_id) const;

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

  /*! Write grid data to a file.
   *
   * \param file_name writes grid file with the prescribed name.
   * \param dir_name writes grid in the prescribed directory.
   * \param show_padding writes grid with or without padding cells.
   */
  void writeToFile(std::string file_name = "grid.dat", std::string dir_name = ".", bool show_padding = false);

 private:
  int m_dimension, m_nx, m_ny, m_nz, m_pad;
  size_t m_total_cells;
  int m_mask[3];
  Vec3<T> m_box_min, m_box_max;
  Vec3<T> m_dx, m_one_over_dx;
  std::vector<Vec3<T>> m_grid;
};

template <typename T, int DIM>
const Mat3<int> Grid<T, DIM>::axis_vectors = Mat3<int>(1, 0, 0, 0, 1, 0, 0, 0, 1);

}  // namespace CPU
}  // namespace GALS
