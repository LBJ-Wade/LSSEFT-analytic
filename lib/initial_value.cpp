//
// Created by David Seery on 22/08/2017.
// --@@
// Copyright (c) 2017 University of Sussex. All rights reserved.
//
// This file is part of the Sussex Effective Field Theory for
// Large-Scale Structure analytic calculation platform (LSSEFT-analytic).
//
// LSSEFT-analytic is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// LSSEFT-analytic is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with LSSEFT-analytic.  If not, see <http://www.gnu.org/licenses/>.
//
// @license: GPL-2
// @contributor: David Seery <D.Seery@sussex.ac.uk>
// --@@
//

#include <sstream>

#include "initial_value.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


initial_value::initial_value(GiNaC::symbol k_, GiNaC::symbol obj_, symbol_factory& sf_)
  : k(std::move(k_)),
    obj(std::move(obj_)),
    sf(sf_)
  {
  }


initial_value_set::initial_value_set(std::initializer_list<initial_value> vals)
  : initial_value_set() // forward to base constructor
  {
    for(const auto& v : vals)
      {
        this->insert(v);
      }
  }


initial_value_set& initial_value_set::insert(initial_value v)
  {
    // check whether the momentum carried by this initial value has already been used
    auto t = this->ks.find(v.get_momentum());
    
    // if momentum already exists then throw an exception
    if(t != this->ks.end())
      {
        std::ostringstream msg;
        msg << ERROR_REPEATED_INITIAL_MOMENTUM << " " << v.get_momentum();
        throw exception(msg.str(), exception_code::initial_value_error);
      }
    
    // insert new momentum into database
    auto r1 = this->ks.insert(v.get_momentum());
    if(!r1.second)
      throw exception(ERROR_INITIAL_VALUE_INSERT_FAILED, exception_code::initial_value_error);
    
    // insert initial value into database
    auto r2 = this->vars.insert(std::move(v));
    if(!r2.second)
      throw exception(ERROR_INITIAL_VALUE_INSERT_FAILED, exception_code::initial_value_error);
    
    return *this;
  }
