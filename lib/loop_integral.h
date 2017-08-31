//
// Created by David Seery on 30/08/2017.
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

#ifndef LSSEFT_ANALYTIC_LOOP_INTEGRAL_H
#define LSSEFT_ANALYTIC_LOOP_INTEGRAL_H


#include <iostream>
#include <services/symbol_factory.h>

#include "utilities/GiNaC_utils.h"


namespace loop_integral_impl
  {

    //! substitution list for Rayleigh momenta
    using subs_list = GiNaC::exmap;

    //! time function
    using time_function = GiNaC::ex;

  }   // namespace loop_integral_impl


//! forward-declare loop_integral
class loop_integral;

//! perform stream insertion
std::ostream& operator<<(std::ostream& str, const loop_integral& obj);


//! loop_integral captures the components needed to construct a component of a 2-point correlation function
class loop_integral
  {

    // TYPES

  public:

    //! pull in time_function
    using time_function = loop_integral_impl::time_function;

    //! pull in subs_list
    using subs_list = loop_integral_impl::subs_list;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    loop_integral(time_function tm_, GiNaC::ex K_, GiNaC::ex ws_, GiNaC_symbol_set lm_, GiNaC_symbol_set em_,
                  subs_list rm_, symbol_factory& sf_)
      : tm(std::move(tm_)),
        K(std::move(K_)),
        WickProduct(std::move(ws_)),
        loop_momenta(std::move(lm_)),
        external_momenta(std::move(em_)),
        Rayleigh_momenta(std::move(rm_)),
        sf(sf_)
      {
      }

    //! destructor is default
    ~loop_integral() = default;


    // TRANSFORMATIONS

  public:

    //! preform reduction of angular integrals
    void reduce_angular_integrals();

  protected:

    //! reduce angular integrals -- one loop implementation
    void reduce_angular_integrals_one_loop();

    //! convert loop momenta to a canonical form
    void canonicalize_loop_labels();

    //! convert Rayleigh momenta to a canonical form
    void canonicalize_Rayleigh_labels();

    //! convert all dot products to cosines
    void dot_products_to_cos();

    //! match arguments of Wick products to Rayleigh momenta
    void match_Wick_to_Rayleigh();
    
    //! convert cosines containing a given vector q to Legendre polynomials
    void cosines_to_Legendre(const GiNaC::symbol& q);

    //! convert Legendre polynomials containing a given vector q to cosines
    void Legendre_to_cosines(const GiNaC::symbol q);


    // SERVICES

  public:

    //! write self to stream
    void write(std::ostream& out) const;


    // INTERNAL DATA

  private:

    // DELEGATES

    //! cache reference to symbol factory
    symbol_factory& sf;


    // KERNEL DATA

    //! time function
    time_function tm;

    //! kernel
    GiNaC::ex K;

    //! string of 2pfs generated by Wick contraction
    GiNaC::ex WickProduct;

    //! set of loop momenta
    GiNaC_symbol_set loop_momenta;

    //! set of external momenta
    GiNaC_symbol_set external_momenta;

    //! set of momenta requiring Rayleigh expansion
    subs_list Rayleigh_momenta;


    // INTEGRATION DATA

    //! set of integration variables
    GiNaC_symbol_set integration_vars;

    //! integration measures
    GiNaC::ex measure{1};


  };


#endif //LSSEFT_ANALYTIC_LOOP_INTEGRAL_H