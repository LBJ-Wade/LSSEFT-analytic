//
// Created by David Seery on 21/08/2017.
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

#ifndef LSSEFT_ANALYTIC_MESSAGES_EN_H
#define LSSEFT_ANALYTIC_MESSAGES_EN_H


constexpr auto WARNING_LABEL = "warning:";

constexpr auto ERROR_SYMBOL_INSERTION_FAILED = "Internal error: symbol insertion failed";
constexpr auto ERROR_INITIAL_VALUE_INSERT_FAILED = "Internal error: initial value insertion failed";

constexpr auto ERROR_KERNEL_NOT_SCALAR = "Kernel is not a scalar function of the momenta";
constexpr auto ERROR_KERNEL_NOT_RATIONAL = "Kernel is not a rational function of the momenta";
constexpr auto ERROR_UNKNOWN_MOMENTA_SING = "Kernel depends on unknown momentum vector";
constexpr auto ERROR_UNKNOWN_MOMENTA_PLURAL = "Kernel depends on unknown momentum vectors";
constexpr auto ERROR_UNUSED_MOMENTA_SING = "Kernel does not depend on available momentum vector";
constexpr auto ERROR_UNUSED_MOMENTA_PLURAL = "Kernel does not depend on available momentum vectors";
constexpr auto ERROR_REPEATED_INITIAL_MOMENTUM = "Attempt to add new initial value with duplicate momentum";

constexpr auto WARNING_KERNEL_EXPRESSION = "Kernel expression";


#endif //LSSEFT_ANALYTIC_MESSAGES_EN_H
