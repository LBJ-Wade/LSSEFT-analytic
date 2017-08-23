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

#ifndef LSSEFT_ANALYTIC_INITIAL_VALUE_H
#define LSSEFT_ANALYTIC_INITIAL_VALUE_H


#include <set>

#include "vector.h"

#include "services/symbol_factory.h"
#include "utilities/GiNaC_utils.h"


//! initial_value represents the initial condition for a stochastic field,
//! eg \delta*_k or similar
class initial_value
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! default constructor is deleted
    initial_value() = delete;
    
    //! destructor is default
    ~initial_value() = default;

  protected:
    
    //! constructor accepts momentum vector and GiNaC symbol specifying
    //! which initial condition this is for; also need a symbol_factory
    //! reference for constructing vectors later
    initial_value(GiNaC::symbol k_, GiNaC::symbol obj_, symbol_factory& _sf);
    
    
    // INTERFACE
    
  public:
    
    //! get momentum vector carried by this variable
    const GiNaC::symbol& get_momentum() const { return this->k; }
    
    //! allow implicit conversion to a vector corresponding to the momentum carried by this field
    operator vector() const { return this->sf.make_vector(this->k); }
    
    //! get GiNaC symbol representing the variable type
    const GiNaC::symbol& get_symbol() const { return this->obj; }
    
    
    // INTERNAL DATA
    
  private:
    
    //! cache momentum vector for this initial value;
    //! note this is a GiNaC symbol rather than a vector instance because we insist that
    //! each stochastic variable has its own single, independent momentum -- they can't
    //! be eg. linear combinations, which would be allowed by a vector
    GiNaC::symbol k;
    
    //! type of object
    GiNaC::symbol obj;
    
    //! capture reference to symbol factory
    symbol_factory& sf;
    
    
    friend class symbol_factory;
  
  };


// specialize std::less to work for initial_value; this allows it to be used in std::set<>
namespace std
  {
    
    // comparison uses only the momenta, not the symbol name
    // this enables us to distinguish fields of equal type, but carrying different momenta
    template<>
    struct less<initial_value>
      {
        bool operator()(const initial_value& a, const initial_value& b) const
          {
            const auto& a_sym = a.get_momentum();
            const auto& b_sym = b.get_momentum();
            
            // have to perform lexical comparison on name ourselves, since GiNaC doesn't
            // automatically provide an ordering on symbols (ie. usually both a < b and b < a,
            // so they compare equal to std::set<>)
            if(a_sym.get_name() < b_sym.get_name()) return true;
            return false;
          }
      };
    
  }   // namespace std


//! initial_value_set represents an unordered set of initial value objects
class initial_value_set
  {
    
    // TYPES
    
  public:
    
    //! database of momentum vectors carried by these initial values
    //! we have to override the normal Compare operator std::less<>,
    //! since GiNaC doesn't provide a sensible ordering on symbols
    //! (any symbol typically compares less than any other symbol, so they
    //! appear equivalent to std::set<>)
    using k_db = GiNaC_symbol_set;
  
  protected:
    
    //! database of initial value objects
    using iv_db = std::set< initial_value >;
    
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor is default
    initial_value_set() = default;
    
    //! initializer_list constructor accepts list of values to insert
    initial_value_set(std::initializer_list<initial_value> vals);
    
    //! destructor is default
    ~initial_value_set() = default;
    
    
    // ITERATORS
    
  public:
    
    // k-database iterators
    k_db::const_iterator          k_cbegin() const      { return this->ks.cbegin(); }
    k_db::const_iterator          k_cend() const        { return this->ks.cend(); }
    k_db::const_reverse_iterator  k_crbegin() const     { return this->ks.crbegin(); }
    k_db::const_reverse_iterator  k_crend() const       { return this->ks.crend(); }
    
    // value database iterators
    iv_db::const_iterator         value_cbegin() const  { return this->vars.cbegin(); }
    iv_db::const_iterator         value_cend() const    { return this->vars.cend(); }
    iv_db::const_reverse_iterator value_crbegin() const { return this->vars.crbegin(); }
    iv_db::const_reverse_iterator value_crend() const   { return this->vars.crend(); }
    
    
    // INTERFACE
    
  public:
    
    //! get the set of momenta used by this set of stochastic variables
    const k_db& get_momenta() const { return this->ks; }
    
    //! add a stochastic variable to the set
    initial_value_set& insert(initial_value v);
    
    
    
    // INTERNAL DATA
    
  private:
    
    //! set of initial value objects
    iv_db vars;
    
    //! set of momenta used by these initial value objects
    k_db ks;
    
  };


#endif //LSSEFT_ANALYTIC_INITIAL_VALUE_H
