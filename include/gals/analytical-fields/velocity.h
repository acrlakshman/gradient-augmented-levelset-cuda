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

#include <map>
#include <string>
#include <vector>

#include "gals/cpu/levelset-velocity.h"
#include "gals/input-fields/velocity.h"
#include "gals/utilities/array.h"
#include "gals/utilities/grid.h"
#include "gals/utilities/vec_n.h"

namespace GALS
{
namespace ANALYTICAL_FIELDS
{
/*! enums for velocity fields.
 */
enum class VelocityFieldNames { NOT_DEFINED, CIRCULAR };

/*! enum maps for velocity fields.
 */
static std::map<std::string, VelocityFieldNames> velocity_name_map{{"CIRCULAR", VelocityFieldNames::CIRCULAR}};

/*! \class Velocity
 *
 * Class to create analytical velocity fields.
 */
template <typename T_GRID, typename T = double>
class Velocity
{
 public:
  /*! Constructor with grid.
   *
   * \param grid grid.
   * \param inputs velocity input fields of type GALS::INPUT_FIELDS::Velocity
   */
  Velocity(const T_GRID& grid, const GALS::INPUT_FIELDS::Velocity& inputs);

  /*! Remove default constructor.
   */
  Velocity() = delete;

  /*! Destructor
   */
  ~Velocity();

  /*! Return grid.
   *
   * \return grid.
   */
  const T_GRID& grid() const { return m_grid; }

  /*! Compute velocity field.
   *
   * \param positions array of positions where velocity needs to be computed.
   * \param time time, for some velocity fields.
   * \param levelset_velocity velocity field which will be updated by this function.
   */
  void compute(const GALS::CPU::Array<T_GRID, GALS::CPU::Vec3<T>>& positions, const T time,
               GALS::CPU::LevelsetVelocity<T_GRID, T>& levelset_velocity);

 private:
  const T_GRID& m_grid;
  const GALS::INPUT_FIELDS::Velocity& m_inputs;
};

}  // namespace ANALYTICAL_FIELDS
}  // namespace GALS
