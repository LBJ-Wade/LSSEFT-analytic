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

#ifndef LSSEFT_ANALYTIC_SYMBOL_FACTORY_H
#define LSSEFT_ANALYTIC_SYMBOL_FACTORY_H


#include <map>
#include <memory>

#include "shared/defaults.h"

#include "boost/optional.hpp"

#include "ginac/ginac.h"


//! forward-declare vector
class vector;


//! symbol factory provides a unified API for building GiNaC symbols
class symbol_factory
  {
    
    // TYPES
    
  protected:
    
    //! type for database key
    using key_type = std::pair< std::string, boost::optional<std::string> >;
    
    //! type for symbol database
    using symbol_db = std::map< key_type, std::unique_ptr<GiNaC::symbol> >;
    
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor accepts default dimension of index-space
    explicit symbol_factory(unsigned int d_ = LSSEFT_DEFAULT_INDEX_DIMENSION);
    
    //! destructor is default
    ~symbol_factory() = default;
    
    
    // OBTAIN SYMBOLS
    
  public:
    
    //! manufacture a GiNaC symbol corresponding to a given name
    GiNaC::symbol& make_symbol(std::string name, boost::optional<std::string> latex_name = boost::none);
    
    //! make a vector object
    vector make_vector(std::string name, boost::optional<std::string> latex_name = boost::none);
    
    
    // INDEX SERVICES
    
  public:
    
    //! make a dummy index
    GiNaC::idx make_dummy_index();
    
    
    // INTERNAL DATA
    
  private:
    
    //! symbol database
    symbol_db symbols;
    
    //! cache dimension of index space
    unsigned int index_dimension;
    
    //! counter for unique index symbols
    unsigned int index_count{0};
  
  };


#endif //LSSEFT_ANALYTIC_SYMBOL_FACTORY_H
