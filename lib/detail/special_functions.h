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

#ifndef LSSEFT_ANALYTIC_SPECIAL_FUNCTIONS_H
#define LSSEFT_ANALYTIC_SPECIAL_FUNCTIONS_H


#include "ginac/ginac.h"


namespace special
  {

    //! spherical bessel function of order nu
    DECLARE_FUNCTION_2P(j);

  }   // namespace special


namespace Angular
  {

    //! cosine between two vectors
    DECLARE_FUNCTION_2P(Cos);

    //! Legendre polynomial of order nu
    DECLARE_FUNCTION_3P(LegP);

  }   // namespace Angular


namespace Fabrikant
  {

    //! declare Fabrikant 3-Bessel integral
    DECLARE_FUNCTION_6P(FabJ)

  }   // namespace Fabrikant


#endif //LSSEFT_ANALYTIC_SPECIAL_FUNCTIONS_H
