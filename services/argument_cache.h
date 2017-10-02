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

#ifndef LSSEFT_ANALYTIC_ARGUMENT_CACHE_H
#define LSSEFT_ANALYTIC_ARGUMENT_CACHE_H


#include "boost/filesystem/operations.hpp"


//! argument cache serves as a central repository for behaviour controls
class argument_cache
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor captures argc, argv parameters from environment
    argument_cache(int argc, char**& argv);

    //! destructor is default
    ~argument_cache() = default;


    // ACCESSORS

  public:

    //! get auto symmetrize status
    bool get_auto_symmetrize() const;

    //! get symmetrize-22 status
    bool get_symmetrize_22() const;

    //! get output root
    const boost::filesystem::path& get_output_root() const;


    // INTERNAL DATA

  private:

    // SYMMETRIZATION

    //! auto-symmetrize kernels?
    bool auto_symmetrize{true};

    //! explicitly symmetrize 22 kernels?
    bool symmetrize_22{true};


    // BACKEND

    //! generate counterterm list?
    bool counterterms{false};
  public:
    bool get_counterterms() const;

  private:

    //! root for output file
    boost::filesystem::path output_root;

  };


#endif //LSSEFT_ANALYTIC_ARGUMENT_CACHE_H
