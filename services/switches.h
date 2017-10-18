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

#ifndef LSSEFT_ANALYTIC_SWITCHES_H
#define LSSEFT_ANALYTIC_SWITCHES_H


constexpr auto SWITCH_VERSION            = "version";
constexpr auto HELP_VERSION              = "show version information";

constexpr auto SWITCH_HELP               = "help";
constexpr auto HELP_HELP                 = "display brief usage information";

constexpr auto SWITCH_COUNTERTERMS       = "counterterms";
constexpr auto SWITCH_NO_COUNTERTERMS    = "no-counterterms";
constexpr auto HELP_COUNTERTERMS         = "compute counterterms";

constexpr auto SWITCH_AUTO_SYMMETRIZE    = "auto-symmetrize";
constexpr auto SWITCH_NO_AUTO_SYMMETRIZE = "no-auto-symmetrize";
constexpr auto HELP_AUTO_SYMMETRIZE      = "automatically symmetrize Fourier kernels";

constexpr auto SWITCH_22_SYMMETRIZE      = "symmetrize-22";
constexpr auto SWITCH_NO_22_SYMMETRIZE   = "no-symmetrize-22";
constexpr auto HELP_22_SYMMETRIZE        = "explicitly symmetrize 22 integrals after angular reduction";

constexpr auto SWITCH_OUTPUT             = "output,o";
constexpr auto SWITCH_OUTPUT_LONG        = "output";
constexpr auto HELP_OUTPUT               = "set output root name";

constexpr auto SWITCH_MATHEMATICA_OUTPUT = "mathematica-output";
constexpr auto HELP_MATHEMATICA_OUTPUT   = "write Mathematica script for loop integrals";


#endif //LSSEFT_ANALYTIC_SWITCHES_H
