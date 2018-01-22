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


constexpr auto FORMAT_HOUR_LABEL = "h";
constexpr auto FORMAT_MINUTE_LABEL = "m";
constexpr auto FORMAT_SECOND_LABEL = "s";

constexpr auto WARNING_LABEL = "warning:";

constexpr auto LABEL_PK = "One-loop power spectrum";
constexpr auto LABEL_PK_TREE = "Tree level";
constexpr auto LABEL_PK_13 = "Loop level 13";
constexpr auto LABEL_PK_22 = "Loop level 22";

constexpr auto ERROR_SYMBOL_INSERTION_FAILED = "Internal error: symbol insertion failed";
constexpr auto ERROR_INITIAL_VALUE_INSERT_FAILED = "Internal error: initial value insertion failed";
constexpr auto ERROR_KERNEL_INSERT_FAILED = "Internal error: kernel insertion failed";
constexpr auto ERROR_KERNEL_COPY_INSERT_FAILED = "Internal error: kernel insertion failed on copy";
constexpr auto ERROR_KERNEL_TRANSFORM_INSERT_FAILED = "Internal error: kernel insertion failed during transformation step";
constexpr auto ERROR_CANT_ADD_PK_INCOMPATIBLE_MOMENTA = "Internal error: can't add power spectra using incompatible momentum variables";
constexpr auto ERROR_LOOP_INTEGRAL_INSERT_FAILED = "Internal error: loop integral insertion failed";
constexpr auto ERROR_LOOP_INTEGRAL_INSERT_AFTER_REDUCTION = "Internal error: attempt to insert loop integral after angular reduction has been performed";
constexpr auto ERROR_ONE_LOOP_ELEMENT_INSERT_FAILED = "Internal error: 1-loop integral element insertion failed";
constexpr auto ERROR_ONE_LOOP_ELEMENT_INSERT_FAILED_RSD = "Internal error: 1-loop integral element insertion failed during RSD filtering";
constexpr auto ERROR_BACKEND_PK_INSERT_FAILED = "Internal error: Pk insertion into backend failed";
constexpr auto ERROR_BACKEND_KERNEL_INSERT_FAILED = "Internal error: kernel insertion into backend failed";

constexpr auto ERROR_KERNEL_NOT_SCALAR = "Kernel is not a scalar function of the momenta";
constexpr auto ERROR_KERNEL_NOT_RATIONAL = "Kernel is not a rational function of the momenta";
constexpr auto ERROR_UNKNOWN_MOMENTA_SING = "Kernel depends on unknown momentum vector";
constexpr auto ERROR_UNKNOWN_MOMENTA_PLURAL = "Kernel depends on unknown momentum vectors";
constexpr auto ERROR_REPEATED_INITIAL_MOMENTUM = "Attempt to add new initial value with duplicate momentum";
constexpr auto ERROR_SUBSTITION_LABEL_NOT_A_SYMBOL = "Kernel substitution list contains a complex element";
constexpr auto ERROR_SUBSTITUTION_LIST_HAS_IV_MOMENTUM = "Kernel substitution list contains initial momentum";
constexpr auto ERROR_MULTIPLY_KERNEL_UNKNOWN_MOMENTUM = "Kernel multiplicand depends on unknown momentum";
constexpr auto ERROR_SUBSTITUTION_RULE_ALREADY_EXISTS = "Redefinition of existing substitution rule";

constexpr auto ERROR_KERNEL_INITIAL_VALUES_DISAGREE = "Internal error: operator += applied to kernels whose initial values disagree";

constexpr auto ERROR_RAYLEIGH_MOMENTUM_IS_ZERO = "Internal error: Rayleigh momentum is zero";

constexpr auto ERROR_EXPONENT_OF_INDEXED_POWER_IS_NOT_NUMERIC = "Exponent in power of indexed object is not numeric";
constexpr auto ERROR_EXPONENT_OF_INDEXED_POWER_IS_NOT_INTEGER = "Exponent in power of indexed object is not an integer";
constexpr auto ERROR_EXPONENT_OF_INDEXED_POWER_IS_NOT_EVEN = "Exponent in power of indexed object is not an even integer";

constexpr auto ERROR_ODD_CONTRACTIONS = "Odd number of fields in contraction: result is zero";
constexpr auto ERROR_UNKNOWN_VERTEX = "Internal error: unknown vertex when traversing graph";
constexpr auto ERROR_CONTRACTION_FAILURE = "Internal error: mismatch when constructing possible initial value contractions";
constexpr auto ERROR_COULD_NOT_ASSIGN_LOOP_MOMENTUM = "Internal error: could not assign loop momentum";
constexpr auto ERROR_COULD_NOT_ASSIGN_EXTERNAL_MOMENTUM = "Internal error: could not assign external momentum";
constexpr auto ERROR_COULD_NOT_EVALUATE_WICK_CONTRACTION = "Internal error: cannot (yet) evaluate a Wick contraction between two unassigned momenta";

constexpr auto ERROR_EXPECTED_ONE_LOOP_RESULT = "Internal error: expected one-loop result";
constexpr auto ERROR_OUTSTANDING_INDICES = "Simplification of some Pk terms was not complete: indices were left over";
constexpr auto ERROR_INCORRECT_SUBMAP_SIZE = "Internal error: unexpected size for substitution map from Wick contractions";

constexpr auto ERROR_CANT_HANDLE_TENSORS = "Cannot yet handle transformations of objects with more than one index";
constexpr auto ERROR_EXPECTED_INDEX_LABEL = "Expected index label";
constexpr auto ERROR_COULD_NOT_FIND_PARTNER_INDEX = "Could not find paired index for dummy summation";
constexpr auto ERROR_EXPECTED_SYMBOL = "Expected dot product factor to be a symbol";
constexpr auto EXPECTED_EXPONENT_TO_BE_NUMERIC = "Expected exponent to be numeric during loop transformation";
constexpr auto EXPECTED_EXPONENT_TO_BE_INTEGER = "Expected exponent to be integer during loop transformation";
constexpr auto EXPECTED_EXPONENT_TO_BE_EVEN_INTEGER = "Expected exponent to be even integer during loop transformation";
constexpr auto ERROR_CANT_MATCH_NEGATIVE_INDEXED_POWER_TO_RAYLEIGH = "Can't match base of negative indexed power to a Rayleigh momentum";
constexpr auto ERROR_INCONSISTENT_INDEXED_SUM = "Inconsistent indices on sum";
constexpr auto ERROR_RAYLEIGH_MOMENTA_POSITIVE_POWER = "Expected Rayleigh momentum to appear with negative power";

constexpr auto ERROR_COULDNT_COLLECT_COS = "Could not collect cosine terms for transformation to Legendre representation";
constexpr auto ERROR_LEGENDRE_TRANSFORM_UNEXPECTED_SIZE = "Internal error: result matrix has unexpected dimensions during conversion to Legendre representation";
constexpr auto ERROR_BADLY_FORMED_WICK_PRODUCT = "Badly formed Wick product";
constexpr auto ERROR_UNKNOWN_WICK_PRODUCT_LABEL = "Wick product contains unknown single momentum label";
constexpr auto ERROR_CANT_MATCH_WICK_TO_RAYLEIGH = "Can't match momentum argument from Wick product to a Rayleigh momentum";
constexpr auto ERROR_BADLY_FORMED_TOP_LEVEL_MOMENTUM_KERNEL = "Badly formed top-level momentum kernel: expected a sum of rational terms";
constexpr auto ERROR_MULTIPLE_RAYLEIGH_MOMENTA_NOT_IMPLEMENTED = "Reduction of multiple Rayleigh momenta is not yet implemented";
constexpr auto ERROR_ONE_LOOP_REDUCE_WITH_TREE = "Internal error: one-loop reduction formula applied to a tree-level term";
constexpr auto ERROR_ONE_LOOP_REDUCE_WITH_MULTIPLE_LOOPS = "Internal error: one-loop reduction formula applied to a multiple-loop term";
constexpr auto ERROR_BADLY_FORMED_TOP_LEVEL_LEGENDRE_SUM = "Badly formed Legendre representations";
constexpr auto ERROR_BADLY_FORMED_LEGENDRE_SUM_TERM = "Badly formed term in Legendre representation";
constexpr auto ERROR_CANT_INTEGRATE_MORE_THAN_TWO_LEGP = "Cannot integrate products of more than two Legendre polynomials";
constexpr auto ERROR_NO_LOOPQ_IN_RAYLEIGH = "Loop momentum does not appear in expression for Rayleigh momentum";
constexpr auto ERROR_DEGREE_OF_LOOPQ_IN_RAYLEIGH_TOO_LARGE = "Loop momentum appears non-linearly in Rayleigh momentum";
constexpr auto ERROR_NO_KEXT_IN_RAYLEIGH = "No external momentum appears in expression for Rayleigh momentum";
constexpr auto ERROR_DEGREE_OF_KEXT_IN_RAYLEIGH_TOO_LARGE = "External momentum appears non-linearly in Rayleigh momentum";
constexpr auto ERROR_TOO_MANY_KEXT_IN_RAYLEIGH = "Too many external momenta appear in Rayleigh momentum";
constexpr auto ERROR_LOOPQ_HAS_WRONG_COEFF_IN_RAYLEIGH = "Loop momentum appears with incorrect coefficient in Rayleigh momentum";
constexpr auto ERROR_KEXT_HAS_WRONG_COEFF_IN_RAYLEIGH = "External momentum appears with incorrect coefficient in Rayleigh momentum";
constexpr auto ERROR_KERNEL_DEPENDS_ON_ANGULAR_RAYLEIGH_MOMENTUM = "Momentum kernel depends on angular part of Rayleigh momentum";
constexpr auto ERROR_TOO_MANY_KEXT_IN_LEGENDRE_SUM = "Too many Legendre polynomials involving loop momentum and external momentum";
constexpr auto ERROR_FAILED_TO_REDUCE_TRIPLE_PRODUCT_LEGP = "Failed to reduce a triple product of Legendre polynomials";
constexpr auto ERROR_CANNOT_ADD_KERNELS_WITH_UNEQUAL_TIME_FUNCTIONS = "Cannot add momentum kernels with differing time functions";

constexpr auto ERROR_FABJ_FIRST_ARG_NONZERO = "First argument of FabJ is not zero";
constexpr auto ERROR_FABJ_SECOND_ARGS_NOT_EQUAL = "Second and third arguments of FabJ are not zero";
constexpr auto ERROR_FABJ_SECOND_ARGS_NEGATIVE = "Second and third arguments of FabJ are negative";

constexpr auto ERROR_INCONSISTENT_KERNEL_FACTOR = "Cannot classify kernel multiplicand";

constexpr auto ERROR_COMPOSE_LOOP_INTEGRAL_MISMATCHING_TYPE = "Internal error: attempt to compose loop integrals of mismatching types";
constexpr auto ERROR_COMPOSE_ONE_LOOP_ELEMENT_MISMATCHING_TYPE = "Internal error: attempt to compose one-loop integral elements of mismatching types";

constexpr auto ERROR_BACKEND_PK_RSD_ALREADY_REGISTERED_A = "A power spectrum with the name";
constexpr auto ERROR_BACKEND_PK_RSD_ALREADY_REGISTERED_B = "is already registered for output";

constexpr auto ERROR_UNKNOWN_GINAC_FUNCTION = "Unknown mathematical function";
constexpr auto ERROR_BACKEND_POW_ARGUMENTS = "Internal error: power function has unexpected number of arguments";
constexpr auto ERROR_BACKEND_PK_ARGUMENTS = "Internal error: Pk correlator has unexpected number of arguments";

constexpr auto WARNING_UNUSED_MOMENTA_SING = "Kernel does not depend on available momentum vector";
constexpr auto WARNING_UNUSED_MOMENTA_PLURAL = "Kernel does not depend on available momentum vectors";
constexpr auto WARNING_ORDER_ZERO_KERNEL = "Ignoring order-zero kernel";
constexpr auto WARNING_KERNEL_EXPRESSION = "Kernel expression";
constexpr auto WARNING_UNEVALUATED_FABRIKANT = "Unevaluated Fabrikant integral of order";
constexpr auto WARNING_PARAMETER_APPEARS_IN_TIME_FUNCTION_RENORMALIZATION = "Encountered normalization of time-dependent function with parameter";
constexpr auto WARNING_PK_RSD_EMPTY = "Empty RSD Pk group for filter pattern";
constexpr auto WARNING_KERNEL_IS_NOT_IR_SAFE = "Detected failure of IR safety for LSSEFT kernel";


#endif //LSSEFT_ANALYTIC_MESSAGES_EN_H
