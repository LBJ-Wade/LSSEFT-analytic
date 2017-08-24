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
#include <functional>

#include "services/symbol_factory.h"
#include "lib/vector.h"
#include "lib/initial_value.h"
#include "lib/fourier_kernel.h"

#include "SPT/time_functions.h"

#include "utilities/GiNaC_utils.h"


int main(int argc, char* argv[])
  {
    symbol_factory sf;
    
    auto k1 = sf.make_vector("k1");
    auto k2 = sf.make_vector("k2");
    
    // generate symbol for magnitude k
    auto k = sf.make_symbol("k");
    
    auto k1_norm = k1.norm_square();
    std::cout << "k1 norm^2 = " << k1_norm << '\n';
    
    auto i = sf.make_unique_index();
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
    
    auto z = sf.get_z();
    GiNaC::ex growth = SPT::D(z);
    std::cout << "growth factor = " << growth << ", derivative dD/dz = " << GiNaC::diff(growth, z) << '\n';
    
    auto deltaq = sf.make_initial_value("delta");
    auto deltar = sf.make_initial_value("delta");
    auto deltas = sf.make_initial_value("delta");
    
    initial_value_set delta1{deltaq};
    initial_value_set delta2{deltaq, deltar};
    initial_value_set delta3{deltaq, deltar, deltas};
    
    vector q = deltaq;
    vector r = deltar;
    vector s = deltas;
    
    GiNaC::ex alpha = dot(q, q+r) / q.norm_square();
    GiNaC::ex beta = dot(q, r) * (q + r).norm_square() / (2 * q.norm_square() * r.norm_square());
    GiNaC::ex gamma = alpha + beta;
    std::cout << "alpha(q,r) = " << alpha << '\n';
    std::cout << "beta(q,r) = " << beta << '\n';
    
    GiNaC::scalar_products sp2;
    GiNaC::symbol qsym = GiNaC::ex_to<GiNaC::symbol>(q.get_expr());
    GiNaC::symbol rsym = GiNaC::ex_to<GiNaC::symbol>(r.get_expr());
    GiNaC::symbol ssym = GiNaC::ex_to<GiNaC::symbol>(s.get_expr());
    sp2.add(qsym, qsym, qsym*qsym);
    sp2.add(rsym, rsym, rsym*rsym);
    sp2.add(ssym, ssym, ssym*ssym);
    std::cout << "alpha simplified = " << simplify_index(alpha, sp2) << '\n';
    
    auto delta = sf.make_fourier_kernel<3>();

    delta.add(SPT::D(z), delta1, 1);
    delta.add(SPT::DA(z), delta2, alpha);
    delta.add(SPT::DB(z), delta2, gamma);
    
    std::cout << delta << '\n';
    std::cout << delta * 2 << '\n';
    std::cout << delta + delta << '\n';
    std::cout << delta - delta << '\n';
    std::cout << delta * delta << '\n';
    
    return EXIT_SUCCESS;
  }
