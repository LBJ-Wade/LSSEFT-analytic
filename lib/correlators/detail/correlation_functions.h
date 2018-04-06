//
// Created by David Seery on 06/04/2018.
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

#ifndef LSSEFT_ANALYTIC_CORRELATION_FUNCTIONS_H
#define LSSEFT_ANALYTIC_CORRELATION_FUNCTIONS_H


#include "ginac/ginac.h"


namespace cfs
  {

    //! declare a GiNaC function to represent a power spectrum
    //! this is really just a placeholder function whose job is to accept three arguments:
    //! the first two specify the correlation function in question and the third specifies
    //! the momentum
    DECLARE_FUNCTION_3P(Pk)

  }   // namespace cfs = correlation functions


#endif //LSSEFT_ANALYTIC_CORRELATION_FUNCTIONS_H
