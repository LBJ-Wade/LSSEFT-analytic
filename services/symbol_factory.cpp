//
// Created by David Seery on 21/08/2017.
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

#include "symbol_factory.h"

#include "lib/vector.h"
#include "lib/initial_value.h"

#include "shared/exceptions.h"
#include "shared/defaults.h"
#include "localizations/messages.h"


symbol_factory::symbol_factory(unsigned int d_)
  : index_dimension(d_),
    z(LSSEFT_REDSHIFT_NAME, LSSEFT_REDSHIFT_LATEX)
  {
  }


const GiNaC::symbol& symbol_factory::make_symbol(std::string name, boost::optional<std::string> latex_name)
  {
    // search for existing definition of this symbol
    auto t = this->symbols.find(std::make_pair(name, latex_name));

    // if a definition was found, return the symbol
    if(t != this->symbols.end()) return t->second;
    
    // if a definition was not found, we need to manufacture a suitable symbol
    GiNaC::symbol sym;
    
    if(latex_name)
      {
        sym = GiNaC::symbol{name, *latex_name};
      }
    else
      {
        sym = GiNaC::symbol{name};
      }
    
    // emplace new element corresponding to this symbol
    key_type key = std::make_pair(std::move(name), std::move(latex_name));
    auto r = this->symbols.emplace(std::move(key), std::move(sym));
    
    // check whether insertion actually occurred
    if(r.second) return r.first->second;
    
    // otherwise, raise an exception
    throw exception(ERROR_SYMBOL_INSERTION_FAILED, exception_code::symbol_error);
  }


const GiNaC::symbol& symbol_factory::get_z()
  {
    return this->z;
  }


GiNaC::idx symbol_factory::make_unique_index()
  {
    // generate unique name for this dummy index
    std::string name = LSSEFT_DEFAULT_INDEX_NAME + std::to_string(this->index_count++);
    
    // generate symbol and convert to an index object
    GiNaC::symbol s{name};
    
    return GiNaC::idx{s, this->index_dimension};
  }


GiNaC::symbol symbol_factory::make_unique_momentum()
  {
    // generate unique name for this momentum variable
    std::string name = LSSEFT_DEFAULT_MOMENTUM_NAME + std::to_string(this->momentum_count++);
    
    // generate symbol
    return GiNaC::symbol{name};
  }


GiNaC::symbol symbol_factory::make_unique_loop_momentum()
  {
    // generate unique name for this loop momentum
    std::string name = LSSEFT_DEFAULT_LOOP_MOMENTUM_NAME + std::to_string(this->loop_count++);

    // generate symbol
    return GiNaC::symbol{name};
  }


vector symbol_factory::make_vector(std::string name, boost::optional<std::string> latex_name)
  {
    auto sym = this->make_symbol(std::move(name), std::move(latex_name));
    
    return this->make_vector(sym);
  }


vector symbol_factory::make_vector(const GiNaC::symbol& s)
  {
    return vector{s, *this};
  }


initial_value symbol_factory::make_initial_value(std::string name, boost::optional<std::string> latex_name)
  {
    auto sym = this->make_symbol(std::move(name), std::move(latex_name));
    auto k = this->make_unique_momentum();
    
    return initial_value{k, sym, *this};
  }


initial_value symbol_factory::make_initial_value(const GiNaC::symbol& s)
  {
    auto k = this->make_unique_momentum();
    
    return initial_value{k, s, *this};
  }
