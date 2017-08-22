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

#include "shared/exceptions.h"
#include "localizations/messages.h"


GiNaC::symbol& symbol_factory::make_symbol(std::string name, boost::optional<std::string> latex_name)
  {
    // search for existing definition of this symbol
    auto t = this->symbols.find(std::make_pair(name, latex_name));

    // if a definition was found, return the symbol
    if(t != this->symbols.end()) return *t->second;
    
    // if a definition was not found, we need to manufacture a suitable symbol
    std::unique_ptr<GiNaC::symbol> sym;
    
    if(latex_name)
      {
        sym = std::make_unique<GiNaC::symbol>(name, *latex_name);
      }
    else
      {
        sym = std::make_unique<GiNaC::symbol>(name);
      }
    
    // emplace new element corresponding to this symbol
    key_type key = std::make_pair(std::move(name), std::move(latex_name));
    auto r = this->symbols.emplace(std::move(key), std::move(sym));
    
    // check whether insertion actually occurred
    if(r.second) return *r.first->second;
    
    // otherwise, raise an exception
    throw exception(ERROR_SYMBOL_INSERTION_FAILED, exception_code::symbol_error);
  }
