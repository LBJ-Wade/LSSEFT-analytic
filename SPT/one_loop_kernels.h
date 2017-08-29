//
// Created by David Seery on 24/08/2017.
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

#ifndef LSSEFT_ANALYTIC_ONE_LOOP_KERNELS_H
#define LSSEFT_ANALYTIC_ONE_LOOP_KERNELS_H


#include "lib/vector.h"


//! alpha kernel
GiNaC::ex alpha(const vector& q, const vector& s, const GiNaC::symbol& eps);

//! beta kernel
GiNaC::ex beta(const vector& q, const vector& s, const GiNaC::symbol& eps);

//! gamma kernel
GiNaC::ex gamma(const vector& q, const vector& s, const GiNaC::symbol& eps);

//! symmetrized alpha kernel
GiNaC::ex alpha_bar(const vector& q, const vector& s, const GiNaC::symbol& eps);

//! symmetrized beta kernel
GiNaC::ex beta_bar(const vector& q, const vector& s, const GiNaC::symbol& eps);

//! symmetrized gamma kernel
GiNaC::ex gamma_bar(const vector& q, const vector& s, const GiNaC::symbol& eps);


#endif //LSSEFT_ANALYTIC_ONE_LOOP_KERNELS_H
