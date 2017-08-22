//
// Created by David Seery on 17/08/2017.
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

#include <iostream>

#include "services/symbol_factory.h"
#include "lib/vector.h"


int main(int argc, char* argv[])
  {
    symbol_factory sf;
    
    auto k1 = sf.make_vector("k1");
    auto k2 = sf.make_vector("k2");
    
    // generate symbol for magnitude k
    auto k = sf.make_symbol("k");
    
    auto k1_norm = k1.norm_square();
    std::cout << "k1 norm^2 = " << k1_norm << '\n';
    
    auto i = sf.make_dummy_index();
    std::cout << "k2 indexed = " << k2[i] << '\n';
    
    auto lin = k1 + 3 * k2;
    std::cout << "lin indexed = " << lin[i] << '\n';
    
    auto dotp = dot(k2, lin);
    std::cout << "k2.lin = " << dotp << '\n';
    
    GiNaC::scalar_products sp;
    sp.add(k1.get_expr(), k1.get_expr(), k);
    sp.add(k2.get_expr(), k2.get_expr(), k);
    sp.add(k1.get_expr(), k2.get_expr(), -k);
    std::cout << "k1 norm^2 simplified = " << k1_norm.simplify_indexed(sp) << '\n';
    std::cout << "k2.lin simplified = " << dotp.simplify_indexed(sp) << '\n';
    
    return EXIT_SUCCESS;
  }
