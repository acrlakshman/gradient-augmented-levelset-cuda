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

#include "gals/input-parser/input-parser-base.h"
#include "gals/input-parser/general.h"
#include "gals/input-parser/grid.h"
#include "gals/input-parser/levelset.h"
#include "gals/input-parser/time.h"
#include "gals/input-parser/velocity.h"

template <typename FIELD>
GALS::INPUT_PARSER::InputParserBase<FIELD>::InputParserBase()
{
}

template <typename FIELD>
GALS::INPUT_PARSER::InputParserBase<FIELD>::~InputParserBase()
{
}

template <typename FIELD>
void GALS::INPUT_PARSER::InputParserBase<FIELD>::parse(const YAML::Node &field,
                                                       GALS::INPUT_FIELDS::InputFields *p_input_fields)
{
  FIELD()(field, p_input_fields);
}

template class GALS::INPUT_PARSER::InputParserBase<GALS::INPUT_PARSER::General>;
template class GALS::INPUT_PARSER::InputParserBase<GALS::INPUT_PARSER::Grid>;
template class GALS::INPUT_PARSER::InputParserBase<GALS::INPUT_PARSER::Time>;
template class GALS::INPUT_PARSER::InputParserBase<GALS::INPUT_PARSER::Velocity>;
template class GALS::INPUT_PARSER::InputParserBase<GALS::INPUT_PARSER::Levelset>;
