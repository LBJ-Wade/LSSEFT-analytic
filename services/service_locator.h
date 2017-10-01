//
// Created by David Seery on 01/10/2017.
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

#ifndef LSSEFT_ANALYTIC_SERVICE_LOCATOR_H
#define LSSEFT_ANALYTIC_SERVICE_LOCATOR_H


#include "argument_cache.h"
#include "symbol_factory.h"


//! forward-declare fourier_kernel
template <unsigned int N>
class fourier_kernel;


//! service_locator is a service locator class
class service_locator
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    service_locator(argument_cache& ac_, symbol_factory& sf_);

    //! destructor is default
    ~service_locator() = default;

    //! disable copying
    service_locator(const service_locator& obj) = delete;


    // FACTORY FUNCTIONS

  public:

    //! make an empty Fourier kernel
    template <unsigned int N>
    fourier_kernel<N> make_fourier_kernel();


    // ACCESSORS

  public:

    //! get argument cache
    argument_cache& get_argunent_cache() { return this->args; }

    //! get symbol factory
    symbol_factory& get_symbol_factory() { return this->sf; }


    // INTERNAL DATA

  public:

    //! capture reference to argument cache
    argument_cache& args;

    //! capture reference to symbol factory
    symbol_factory& sf;

  };


template <unsigned int N>
fourier_kernel<N> service_locator::make_fourier_kernel()
  {
    return fourier_kernel<N>{*this};
  }


#endif //LSSEFT_ANALYTIC_SERVICE_LOCATOR_H
