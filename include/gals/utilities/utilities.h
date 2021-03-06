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

#include <math.h>
#include <string.h>

namespace GALS
{
#ifndef BUILD_TESTS
#define GALS_FUNCTION_NOT_IMPLEMENTED(var)                        \
  std::cout << "FUNCTION_NOT_IMPLEMENTED: " << #var << std::endl; \
  exit(0);
#else
#define GALS_FUNCTION_NOT_IMPLEMENTED(var) std::cout << "FUNCTION_NOT_IMPLEMENTED: " << #var << std::endl;
#endif

#ifndef BUILD_TESTS
#define GALS_ABORT(var)                        \
  std::cout << "ERROR: " << #var << std::endl; \
  exit(0);
#else
#define GALS_ABORT(var) std::cout << "ERROR: " << #var << std::endl;
#endif

static double pi() { return atan(1) * 4; }

static double VSMALL = 1e-10;
static const double one_third = 1. / 3.;
static const double two_thirds = 2. / 3.;

/*! Check for equality within a tolerance.
 *
 * Tolerance is by default 1e-10.
 *
 * \param a first argument
 * \param b second argument
 *
 * \return number of elements.
 */
static bool is_equal(double a, double b) { return fabs(a - b) <= VSMALL; }

//! Check for equality between integers.
static bool is_equal(int a, int b) { return a == b; }

//! Compute square.
template <typename T>
static T sqr(T a)
{
  return a * a;
}

//! Compute cube.
template <typename T>
static T cube(T a)
{
  return a * a * a;
}

}  // namespace GALS
