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

#ifndef LSSEFT_ANALYTIC_DEFAULTS_H
#define LSSEFT_ANALYTIC_DEFAULTS_H


//! default dimension of spatial indices = 3
constexpr unsigned int LSSEFT_DEFAULT_INDEX_DIMENSION = 3;

//! default name for dummy indices
constexpr auto LSSEFT_DEFAULT_INDEX_NAME = "i";

//! default name for dummy momentum
constexpr auto LSSEFT_DEFAULT_MOMENTUM_NAME = "q";

//! default name for dummy loop momentum
constexpr auto LSSEFT_DEFAULT_LOOP_MOMENTUM_NAME = "l";

//! default name for Rayliegh momenta
constexpr auto LSSEFT_DEFAULT_RAYLEIGH_MOMENTUM_NAME = "r";

//! default text name for redshift variable
constexpr auto LSSEFT_REDSHIFT_NAME = "z";

//! default LaTeX name for redshift variable
constexpr auto LSSEFT_REDSHIFT_LATEX = "z";

//! default text name for regulator (needed for denominators that can become zero)
constexpr auto LSSEFT_EPSILON_NAME = "eps";

//! default LaTeX name for regulator (needed for denominators that can become zero)
constexpr auto LSSEFT_EPSILON_LATEX = "\epsilon";


#endif //LSSEFT_ANALYTIC_DEFAULTS_H
