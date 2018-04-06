//
// Created by David Seery on 29/08/2017.
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

#include "relabel_product.h"


GiNaC::ex relabel_index_product(const GiNaC::ex& a, const GiNaC::ex& b, service_locator& loc)
  {
    GiNaC::exmap idx_map;
    auto& sf =  loc.get_symbol_factory();

    // need only relabel indices that occur >= 2 times; can allow single occurrences to be contracted
    const auto& a_idxs = get_expr_indices(a, 2);
    const auto& b_idxs = get_expr_indices(b, 2);

    for(const auto& idx : b_idxs)
    {
      // does this index already exist in the a set? if not, nothing to do
      if(a_idxs.find(idx) == a_idxs.end()) continue;

      // otherwise, manufacture a new index
      // note that we remap the index *value* (ie. its symbolic label), not the index itself
        const auto relabel = sf.make_unique_index();
        idx_map[idx] = GiNaC::ex_to<GiNaC::symbol>(relabel.get_value());
    }

    return a * b.subs(idx_map);
  }
