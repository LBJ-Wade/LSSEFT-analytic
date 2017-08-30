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
#include "lib/initial_value.h"
#include "lib/fourier_kernel.h"


//! alpha kernel
fourier_kernel_impl::kernel alpha(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf);

//! beta kernel
fourier_kernel_impl::kernel beta(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf);

//! gamma kernel
fourier_kernel_impl::kernel gamma(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf);

//! symmetrized alpha kernel
fourier_kernel_impl::kernel alpha_bar(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf);

//! symmetrized beta kernel
fourier_kernel_impl::kernel beta_bar(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf);

//! symmetrized gamma kernel
fourier_kernel_impl::kernel gamma_bar(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf);


#endif //LSSEFT_ANALYTIC_ONE_LOOP_KERNELS_H
