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
#include "lib/Pk_one_loop.h"

#include "SPT/time_functions.h"
#include "SPT/one_loop_kernels.h"



int main(int argc, char* argv[])
  {
    symbol_factory sf;
    
    // redshift z is the time variable
    auto z = sf.get_z();
    
    // r is the unit line-of-sight vector to Earth
    auto r_sym = sf.make_symbol("r");
    auto r = sf.make_vector(r_sym);
    
    // mu is RSD parameter = r.\hat{k} = r.k / |k|
    auto mu = sf.make_symbol("mu");
    
    // manufacture placeholder stochastic initial values delta*_q, delta*_s, delta*_t
    // (recall we skip delta*_r because r is also the line-of-sight variable)
    auto deltaq = sf.make_initial_value("delta");
    auto deltas = sf.make_initial_value("delta");
    auto deltat = sf.make_initial_value("delta");
    
    // linear set is delta*_q
    initial_value_set iv_q{deltaq};
    
    // quadratic set is delta*_q delta*_s or delta*_s delta*_t
    initial_value_set iv_qs{deltaq, deltas};
    initial_value_set iv_st{deltas, deltat};
    
    // cubic set is delta*_q delta*_s delta*_t
    initial_value_set iv_qst{deltaq, deltas, deltat};

    // extract momentum vectors from these initial value placeholders
    vector q = deltaq;
    vector s = deltas;
    vector t = deltat;
    
    
    // set up kernels for the dark matter overdensity \delta
    auto delta = sf.make_fourier_kernel<3>();

    // linear order
    delta.add(SPT::D(z), iv_q, 1);
    
    // second order
    delta.add(SPT::DA(z) * alpha(q, s, kernel{iv_qs, sf}, sf));
    delta.add(SPT::DB(z) * gamma(q, s, kernel{iv_qs, sf}, sf));

    // quadratic order
    delta.add((SPT::DD(z) - SPT::DJ(z)) * 2*gamma_bar(s+t, q, alpha_bar(s, t, kernel{iv_qst, sf}, sf), sf));
    delta.add(SPT::DE(z)                * 2*gamma_bar(s+t, q, gamma_bar(s, t, kernel{iv_qst, sf}, sf), sf));
    delta.add((SPT::DF(z) + SPT::DJ(z)) * 2*alpha_bar(s+t, q, alpha_bar(s, t, kernel{iv_qst, sf}, sf), sf));
    delta.add(SPT::DG(z)                * 2*alpha_bar(s+t, q, gamma_bar(s, t, kernel{iv_qst, sf}, sf), sf));
    delta.add(SPT::DJ(z)                * (alpha(s+t, q, gamma_bar(s, t, kernel{iv_qst, sf}, sf), sf)
                                           - 2*alpha(s+t, q, alpha_bar(s, t, kernel{iv_qst, sf}, sf), sf)));
    
    
    // compute kernels for the velocity potential \phi, v = grad phi -> v(k) = i k phi
    auto delta1 = delta.order(1);
    auto delta2 = delta.order(2);
    auto delta3 = delta.order(3);
    
    auto phi1 = InverseLaplacian(-diff_t(delta1));
    auto phi2 = InverseLaplacian(-diff_t(delta2) - delta1*Laplacian(phi1) - gradgrad(phi1, delta1));
    auto phi3 = InverseLaplacian(-diff_t(delta3)
                                 - delta1*Laplacian(phi2) - delta2*Laplacian(phi1)
                                 - gradgrad(phi1, delta2) - gradgrad(phi2, delta1));
    
    auto phi = phi1 + phi2 + phi3;
    
    
    // set up momentum label k
    auto k = sf.make_symbol("k");

    sf.declare_parameter(r_sym);
    sf.declare_parameter(k);
    sf.declare_parameter(mu);
    
    // build expression for the redshift-space \delta
    GiNaC::ex H = FRW::Hub(z);
    auto k1mu = k*mu;
    auto k2mu = -k*mu;
    
    auto r_dot_v = dotgrad(r, phi);

    auto make_delta_rsd = [&](const decltype(k1mu)& kmu) -> auto
      {
        return delta
               - (GiNaC::I / H) * kmu * r_dot_v
               - (GiNaC::I / H) * kmu * (r_dot_v * delta)
               + (GiNaC::ex(1) / (2*H*H)) * kmu*kmu * (r_dot_v * r_dot_v)
               + (GiNaC::ex(1) / (2*H*H)) * kmu*kmu * (r_dot_v * r_dot_v * delta)
               + (GiNaC::I / (3*H*H*H)) * kmu*kmu*kmu * (r_dot_v * r_dot_v * r_dot_v);
      };
    auto delta_rsd_k1 = make_delta_rsd(k1mu);
    auto delta_rsd_k2 = make_delta_rsd(k2mu);
    
    // construct 1-loop \delta power spectrum
    Pk_one_loop Pk_delta{delta, delta, k, sf};
    
    const auto& tree = Pk_delta.get_tree();
    std::cout << "Tree-level P(k):" << '\n';
    std::cout << tree << '\n';

    const auto& P13 = Pk_delta.get_13();
    std::cout << "Loop-level 13 P(k):" << '\n';
    std::cout << P13 << '\n';

    const auto& P22 = Pk_delta.get_22();
    std::cout << "Loop-level 22 P(k):" << '\n';
    std::cout << P22 << '\n';
    
    return EXIT_SUCCESS;
  }
