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


fourier_kernel_impl::kernel alpha(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf)
  {
    fourier_kernel_impl::subs_list Rayleigh_list;
    
    auto R_sym = sf.make_unique_Rayleigh_momentum();
    auto R = sf.make_vector(R_sym);
    
    GiNaC::ex K = dot(q, q+s) / R.norm_square();
    Rayleigh_list[R_sym] = q.get_expr();
    
    return fourier_kernel_impl::kernel{K, iv, GiNaC::ex(1), Rayleigh_list, sf};
  }


fourier_kernel_impl::kernel beta(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf)
  {
    fourier_kernel_impl::subs_list Rayleigh_list;
    
    auto R_sym = sf.make_unique_Rayleigh_momentum();
    auto R = sf.make_vector(R_sym);
    
    auto S_sym = sf.make_unique_Rayleigh_momentum();
    auto S = sf.make_vector(S_sym);

    GiNaC::ex K = dot(q, s) * (q+s).norm_square() / (2 * R.norm_square() * S.norm_square());
    Rayleigh_list[R_sym] = q.get_expr();
    Rayleigh_list[S_sym] = s.get_expr();
    
    return fourier_kernel_impl::kernel{K, iv, GiNaC::ex(1), Rayleigh_list, sf};
  }


fourier_kernel_impl::kernel gamma(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf)
  {
    return alpha(q, s, iv, sf) + beta(q, s, iv, sf);
  }


fourier_kernel_impl::kernel alpha_bar(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf)
  {
    return (alpha(q, s, iv, sf) + alpha(s, q, iv, sf)) / 2;
  }


fourier_kernel_impl::kernel beta_bar(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf)
  {
    return (beta(q, s, iv, sf) + beta(s, q, iv, sf)) / 2;
  }


fourier_kernel_impl::kernel gamma_bar(const vector& q, const vector& s, const initial_value_set& iv, symbol_factory& sf)
  {
    return (gamma(q, s, iv, sf) + gamma(s, q, iv, sf)) / 2;
  }
