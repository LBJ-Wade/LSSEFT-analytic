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

#include "shared/defaults.h"
#include "utilities/GiNaC_utils.h"

#include "boost/optional.hpp"

#include "ginac/ginac.h"


//! forward-declare vector
class vector;

//! forward-declare initial_value
class initial_value;


//! symbol factory provides a unified API for building GiNaC symbols
class symbol_factory
  {
    
    // TYPES
    
  protected:
    
    //! type for database key
    using key_type = std::pair< std::string, boost::optional<std::string> >;
    
    //! type for symbol database
    using symbol_db = std::map< key_type, GiNaC::symbol >;
    
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor accepts default dimension of index-space
    explicit symbol_factory(unsigned int d_ = LSSEFT_DEFAULT_INDEX_DIMENSION);
    
    //! destructor is default
    ~symbol_factory() = default;
    
    
    // OBTAIN SYMBOLS
    
  public:
    
    //! get GiNaC symbol representing redshift z
    const GiNaC::symbol& get_z();

    //! manufacture a GiNaC symbol corresponding to a given name
    const GiNaC::symbol& make_symbol(std::string name, boost::optional<std::string> latex_name = boost::none);
    
    
    // FACTORY FUNCTIONS
    
  public:
    
    //! make a vector object from a name
    vector make_vector(std::string name, boost::optional<std::string> latex_name = boost::none);
    
    //! make a vector object from an existing symbol
    vector make_vector(const GiNaC::symbol& s);
    
    //! make an initial value object
    initial_value make_initial_value(std::string name, boost::optional<std::string> latex_name = boost::none);
    
    //! make an initial value object from an existing symbol
    initial_value make_initial_value(const GiNaC::symbol& s);

    
    // SERVICES
    
  public:
    
    //! make a unique index
    GiNaC::idx make_unique_index();
    
    //! make a unique momentum variable
    GiNaC::symbol make_unique_momentum();

    //! make a unique loop momentum variable
    GiNaC::symbol make_unique_loop_momentum();
    
    //! make a unique Rayleigh variable
    GiNaC::symbol make_unique_Rayleigh_momentum();

    //! make a canonicalized loop momentum
    GiNaC::symbol make_canonical_loop_momentum(unsigned int count);

    //! make a canonicalized Rayleigh momentum
    GiNaC::symbol make_canonical_Rayleigh_momentum(unsigned int count);

    // PARAMETER HANDLING

  public:

    //! declare a symbol to be a parameter
    void declare_parameter(const GiNaC::symbol& s);

    //! get parameter set
    const GiNaC_symbol_set& get_parameters() const { return this->parameters; }
    
    
    // INTERNAL DATA
    
  private:

    // DATABASES

    //! symbol database
    symbol_db symbols;

    //! parameter database; symbols declared as parameters can be included in momentum kernels
    GiNaC_symbol_set parameters;


    // INTERNAL STATE
    
    //! cache dimension of index space
    unsigned int index_dimension;
    
    //! counter for unique index symbols
    unsigned int index_count{0};
    
    //! counter for unique momentum symbols
    unsigned int momentum_count{0};

    //! counter for unique loop momentum symbols
    unsigned int loop_count{0};
    
    //! counter for unique Rayleigh momentum symbols
    unsigned int Rayleigh_count{0};
    
    
    // RESERVED SYMBOLS
    
    //! redshift z
    GiNaC::symbol z;

  };


#endif //LSSEFT_ANALYTIC_SYMBOL_FACTORY_H
