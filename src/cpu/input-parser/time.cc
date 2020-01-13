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

#include "gals/input-parser/time.h"
#include "gals/input-fields/time.h"

GALS::INPUT_PARSER::Time::Time() {}

GALS::INPUT_PARSER::Time::~Time() {}

void GALS::INPUT_PARSER::Time::parse(const YAML::Node &field, GALS::INPUT_FIELDS::InputFields *p_input_fields)
{
  auto &input_fields = *p_input_fields;

  // Parse time section.
  input_fields.m_time->start = field["start"].as<double>();
  input_fields.m_time->end = field["end"].as<double>();
  input_fields.m_time->dt = field["dt"].as<double>();
  input_fields.m_time->constant_dt = field["constant_dt"].as<std::string>();
  input_fields.m_time->cfl_max = field["cfl_max"].as<double>();
  input_fields.m_time->write_interval = field["write_interval"].as<double>();
}

void GALS::INPUT_PARSER::Time::operator()(const YAML::Node &field, GALS::INPUT_FIELDS::InputFields *p_input_fields)
{
  this->parse(field, p_input_fields);
}
