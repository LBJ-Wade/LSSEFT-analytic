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

#include "one_loop_kernels.h"


kernel alpha(const vector& q, const vector& s, kernel ker, symbol_factory& sf)
  {
    auto Q_sym = sf.make_unique_Rayleigh_momentum();
    auto Q = sf.make_vector(Q_sym);
    
    GiNaC::ex K = dot(q, q+s) / Q.norm_square();

    ker.multiply_kernel(K, Q_sym, q.get_expr());
    return ker;
  }


kernel beta(const vector& q, const vector& s, kernel ker, symbol_factory& sf)
  {
    auto Q_sym = sf.make_unique_Rayleigh_momentum();
    auto Q = sf.make_vector(Q_sym);
    
    auto S_sym = sf.make_unique_Rayleigh_momentum();
    auto S = sf.make_vector(S_sym);

    GiNaC::ex K = dot(q, s) * (q+s).norm_square() / 2;
    ker.multiply_kernel(K);
    ker.multiply_kernel(GiNaC::ex(1) / Q.norm_square(), Q_sym, q.get_expr());
    ker.multiply_kernel(GiNaC::ex(1) / S.norm_square(), S_sym, s.get_expr());
    return ker;
  }


kernel gamma(const vector& q, const vector& s, kernel ker, symbol_factory& sf)
  {
    kernel ker_copy = ker;
    return alpha(q, s, ker, sf) + beta(q, s, ker_copy, sf);
  }


kernel alpha_bar(const vector& q, const vector& s, kernel ker, symbol_factory& sf)
  {
    kernel ker_copy = ker;
    return (alpha(q, s, ker, sf) + alpha(s, q, ker_copy, sf)) / 2;
  }


kernel beta_bar(const vector& q, const vector& s, kernel ker, symbol_factory& sf)
  {
    kernel ker_copy = ker;
    return (beta(q, s, ker, sf) + beta(s, q, ker_copy, sf)) / 2;
  }


kernel gamma_bar(const vector& q, const vector& s, kernel ker, symbol_factory& sf)
  {
    kernel ker_copy = ker;
    return (gamma(q, s, ker, sf) + gamma(s, q, ker_copy, sf)) / 2;
  }
