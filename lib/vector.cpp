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
#include "vector.h"


vector::vector(GiNaC::symbol& k, symbol_factory& sf_)
  : expr(k),
    sf(sf_)
  {
  }


vector::vector(GiNaC::ex e, symbol_factory& sf_)
  : expr(std::move(e)),
    sf(sf_)
  {
  }


vector operator-(const vector& a)
  {
    return vector{-a.expr, a.sf};
  }


vector operator-(const vector& a, const vector& b)
  {
    return vector{a.expr - b.expr, a.sf};
  }


vector operator+(const vector& a, const vector& b)
  {
    return vector{a.expr + b.expr, a.sf};
  }


vector operator/(const vector& a, const GiNaC::numeric b)
  {
    return vector{a.expr / b, a.sf};
  }


vector operator*(const vector& a, const GiNaC::numeric b)
  {
    return vector{a.expr * b, a.sf};
  }


vector operator*(const GiNaC::numeric a, const vector& b)
  {
    return vector{a * b.expr, b.sf};
  }


GiNaC::ex dot(const vector a, const vector b)
  {
    auto dummy_index = a.sf.make_unique_index();
    
    auto a_idx = a[dummy_index];
    auto b_idx = b[dummy_index];
    
    GiNaC::ex ab = a_idx * b_idx;
    
    return ab.expand(GiNaC::expand_options::expand_indexed);
  }


GiNaC::ex vector::operator[](const GiNaC::idx& i) const
  {
    GiNaC::ex cmp = GiNaC::indexed(this->expr, i);
    return cmp.expand(GiNaC::expand_options::expand_indexed).simplify_indexed();
  }


GiNaC::ex vector::norm_square() const
  {
    auto dummy_index = this->sf.make_unique_index();
    
    auto expr_idx = this->operator[](dummy_index);
    GiNaC::ex expr_idx_sq = expr_idx * expr_idx;
    
    return expr_idx_sq.expand(GiNaC::expand_options::expand_indexed).simplify_indexed();
  }


GiNaC::ex vector::norm() const
  {
    return GiNaC::sqrt(this->norm_square());
  }
