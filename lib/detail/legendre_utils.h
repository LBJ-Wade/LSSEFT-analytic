//
// Created by David Seery on 31/08/2017.
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

#ifndef LSSEFT_ANALYTIC_LEGENDRE_UTILS_H
#define LSSEFT_ANALYTIC_LEGENDRE_UTILS_H


#include "utilities/GiNaC_utils.h"


//! calculate a given Legendre polynomial
GiNaC::ex LegP(unsigned int n, const GiNaC::ex& t);

//! convert all dot products to cosines
GiNaC::ex dot_products_to_cos(const GiNaC::ex& expr);

//! convert cosines containing a given vector q to Legendre polynomials
GiNaC::ex cosines_to_Legendre(const GiNaC::ex& expr, const GiNaC::symbol& q);

//! convert Legendre polynomials containing a given vector q to cosines
GiNaC::ex Legendre_to_cosines(GiNaC::ex expr, const GiNaC::symbol q);

//! get set of symbols appearing in cosine functions with q
GiNaC_symbol_set get_Cos_pairs(const GiNaC::symbol& q, const GiNaC::ex& expr);

//! get set of symbols appears in LegP functions with q
GiNaC_symbol_set get_LegP_pairs(const GiNaC::ex& expr, const GiNaC::symbol& q);


#endif //LSSEFT_ANALYTIC_LEGENDRE_UTILS_H
