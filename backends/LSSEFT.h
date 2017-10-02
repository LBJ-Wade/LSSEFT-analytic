//
// Created by David Seery on 02/10/2017.
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
#ifndef LSSEFT_ANALYTIC_LSSEFT_H
#define LSSEFT_ANALYTIC_LSSEFT_H


#include <map>
#include <functional>

#include "lib/Pk_rsd.h"


//! LSSEFT is a backend class capable of writing out C++ to implement
//! a given correlation function in the Sussex LSSEFT platform
class LSSEFT
  {

    // TYPES

  protected:

    //! main database
    using db_type = std::map< std::string, std::reference_wrapper<const Pk_rsd> >;

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    LSSEFT() = default;

    //! destructor
    ~LSSEFT() = default;


    // INTERFACE

  public:

    //! add a new power spectrum
    LSSEFT& add(const Pk_rsd& P, std::string name);


    // INTERNAL DATA

  private:

    //! main database
    db_type db;

  };


#endif //LSSEFT_ANALYTIC_LSSEFT_H
