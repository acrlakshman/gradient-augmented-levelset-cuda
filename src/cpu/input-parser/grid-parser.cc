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

#include "grid-parser.h"
#include "../input-fields/grid.h"

#include <vector>

GALS::CPU::GridParser::GridParser() {}

GALS::CPU::GridParser::~GridParser() {}

void GALS::CPU::GridParser::parse(const YAML::Node &field, GALS::INPUT_FIELDS::InputFields *p_input_fields)
{
  auto &input_fields = *p_input_fields;

  // Domain bounds.
  input_fields.m_grid->x_min = field["box"]["x_min"].as<double>();
  input_fields.m_grid->x_max = field["box"]["x_max"].as<double>();
  input_fields.m_grid->y_min = field["box"]["y_min"].as<double>();
  input_fields.m_grid->y_max = field["box"]["y_max"].as<double>();
  input_fields.m_grid->z_min = field["box"]["z_min"].as<double>();
  input_fields.m_grid->z_max = field["box"]["z_max"].as<double>();

  // Number of cells.
  input_fields.m_grid->nx = field["cells"].as<std::vector<int>>()[0];
  input_fields.m_grid->ny = field["cells"].as<std::vector<int>>()[1];
  input_fields.m_grid->nz = field["cells"].as<std::vector<int>>()[2];
}

void GALS::CPU::GridParser::operator()(const YAML::Node &field, GALS::INPUT_FIELDS::InputFields *p_input_fields)
{
  this->parse(field, p_input_fields);
}
