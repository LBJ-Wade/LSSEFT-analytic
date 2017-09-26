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
#include "lib/Pk_rsd.h"
#include "lib/detail/special_functions.h"

#include "SPT/time_functions.h"
#include "SPT/one_loop_kernels.h"


std::vector<std::string> generate_UV_limit(const Pk_rsd_group& group, const GiNaC::symbol& k, unsigned int max_mu, unsigned int max_k)
  {
    std::vector<std::string> output;

    const auto UV_limit = group.get_UV_limit(2*max_k);

    for(unsigned int i = 0; i <= max_k; ++i)
      {
        const auto kval = 2*i;

        std::ostringstream header;
        header << "@ k^" << kval;
        output.push_back(header.str());

        for(unsigned int j = 0; j <= max_mu; ++j)
          {
            const auto muval = 2*j;

            std::ostringstream msg;
            auto expr = GiNaC::collect_common_factors(UV_limit[j].expand().coeff(k, kval));
            msg << "   -> mu^" << muval << " = " << expr;
            output.push_back(msg.str());
          }
      }

    return output;
  }


std::vector<std::string> generate_map(const Pk_rsd_group& group, const GiNaC::symbol& k, unsigned int max_mu, unsigned int max_k)
  {
    std::vector<std::string> output;

    const auto UV_limit = group.get_UV_limit(2*max_k);
    const auto time_funcs = group.get_time_functions();

    for(unsigned int i = 0; i <= max_mu; ++i)
      {
        const auto muval = 2*i;

        std::ostringstream msg;
        unsigned int count = 0;

        for(unsigned int j = 0; j <= max_k; ++j)
          {
            const auto power = 2*j;
            auto expr = UV_limit[i].expand().coeff(k, power);
            if(expr != 0)
              {
                ++count;
                msg << " k^" << power;
              }
          }

        if(count > 0)
          {
            std::ostringstream line;

            line << "   -- mu^" << muval << " ->" << msg.str();
            if(!time_funcs[i].empty())
              {
                line << " (" << time_funcs[i].size() << " time function" << (time_funcs[i].size() != 1 ? "s" : "") << ")";
              }
            output.push_back(line.str());

            if(!time_funcs[i].empty())
              {
                for(unsigned int m = 0; m < time_funcs[i].size(); ++m)
                  {
                    std::ostringstream tline;
                    tline << "   -- -- " << m << ". " << time_funcs[i][m] << '\n';
                    output.push_back(tline.str());
                  }
              }
          }
      }

    return output;
  }


std::vector<std::string> mixing_map(const Pk_rsd& rsd, const GiNaC::symbol& k)
  {
    const auto& rsd_13 = rsd.get_13();
    return generate_map(rsd_13, k, 4, 1);
  }


std::vector<std::string> stochastic_map(const Pk_rsd& rsd, const GiNaC::symbol& k)
  {
    const auto& rsd_22 = rsd.get_22();
    return generate_map(rsd_22, k, 4, 2);
  }


std::vector<std::string> mixing_divergences(const Pk_rsd& rsd, const GiNaC::symbol& k)
  {
    const auto& rsd_13 = rsd.get_13();
    return generate_UV_limit(rsd_13, k, 4, 1);
  }


std::vector<std::string> stochastic_divergences(const Pk_rsd& rsd, const GiNaC::symbol& k)
  {
    const auto& rsd_22 = rsd.get_22();
    return generate_UV_limit(rsd_22, k, 4, 2);
  }


template <typename MapGenerator>
void write_map(const Pk_rsd& Pk_nobias, const Pk_rsd& Pk_b1, const Pk_rsd& Pk_b2, const Pk_rsd& Pk_b3,
               const Pk_rsd& Pk_bG2, const Pk_rsd& Pk_bdG2, const Pk_rsd& Pk_bGamma3,
               const Pk_rsd& Pk_b1b1, const Pk_rsd& Pk_b1b2, const Pk_rsd& Pk_b1b3, const Pk_rsd& Pk_b2b2,
               const Pk_rsd& Pk_b1bG2, const Pk_rsd& Pk_bG2bG2, const Pk_rsd& Pk_b2bG2, const Pk_rsd& Pk_b1bdG2,
               const Pk_rsd& Pk_b1bGamma3, const GiNaC::symbol& k, MapGenerator make_map)
  {
    auto write = [&](std::string name, const Pk_rsd& rsd) -> void
      {
        const auto output = make_map(rsd, k);
        if(!output.empty())
          {
            std::cout << "-- " << name << '\n';
            for(const auto& line : output)
              {
                std::cout << line << '\n';
              }
            std::cout << '\n';
          }
      };

    write("no bias", Pk_nobias);

    write("b1", Pk_b1);
    write("b2", Pk_b2);
    write("b3", Pk_b3);
    write("bG2", Pk_bG2);
    // bG3 gives zero
    write("bdG2", Pk_bdG2);
    write("bGamma3", Pk_bGamma3);

    write("b1 x b1", Pk_b1b1);
    write("b1 x b2", Pk_b1b2);
    write("b1 x b3", Pk_b1b3);
    write("b2 x b2", Pk_b2b2);

    write("b1 x bG2", Pk_b1bG2);
    write("b2 x bG2", Pk_b2bG2);
    write("bG2 x bG2", Pk_bG2bG2);
    // b1 x bG3 gives zero
    write("b1 x bdG2", Pk_b1bdG2);
    write("b1 x bGamma3", Pk_b1bGamma3);
  }


void count_kernels(const Pk_rsd& Pk_nobias, const Pk_rsd& Pk_b1, const Pk_rsd& Pk_b2, const Pk_rsd& Pk_b3,
                   const Pk_rsd& Pk_bG2, const Pk_rsd& Pk_bdG2, const Pk_rsd& Pk_bGamma3,
                   const Pk_rsd& Pk_b1b1, const Pk_rsd& Pk_b1b2, const Pk_rsd& Pk_b1b3, const Pk_rsd& Pk_b2b2,
                   const Pk_rsd& Pk_b1bG2, const Pk_rsd& Pk_bG2bG2, const Pk_rsd& Pk_b2bG2, const Pk_rsd& Pk_b1bdG2,
                   const Pk_rsd& Pk_b1bGamma3)
  {
    size_t kernels = 0;

    auto count = [&](std::string name, const Pk_rsd& group) -> void
      {
        const auto& P13 = group.get_13();
        const auto& P22 = group.get_22();

        const auto P13_time = P13.get_time_functions();
        const auto P22_time = P22.get_time_functions();

        size_t kernels_13 = 0;
        for(unsigned int i = 0; i < P13_time.size(); ++i)
          {
            kernels_13 += P13_time[i].size();
          }

        size_t kernels_22 = 0;
        for(unsigned int i = 0; i < P22_time.size(); ++i)
          {
            kernels_22 += P22_time[i].size();
          }

        std::cout << "-- " << name << " -> 13 kernels = " << kernels_13 << ", 22 kernels = " << kernels_22 << '\n';

        kernels += kernels_13 + kernels_22;
      };

    count("no bias", Pk_nobias);

    count("b1", Pk_b1);
    count("b2", Pk_b2);
    count("b3", Pk_b3);
    count("bG2", Pk_bG2);
    // bG3 gives zero
    count("bdG2", Pk_bdG2);
    count("bGamma3", Pk_bGamma3);

    count("b1 x b1", Pk_b1b1);
    count("b1 x b2", Pk_b1b2);
    count("b1 x b3", Pk_b1b3);
    count("b2 x b2", Pk_b2b2);

    count("b1 x bG2", Pk_b1bG2);
    count("b2 x bG2", Pk_b2bG2);
    count("bG2 x bG2", Pk_bG2bG2);
    // b1 x bG3 gives zero
    count("b1 x bdG2", Pk_b1bdG2);
    count("b1 x bGamma3", Pk_b1bGamma3);

    std::cout << '\n' << "TOTAL KERNELS = " << kernels << '\n';
  }


int main(int argc, char* argv[])
  {
    symbol_factory sf;

    // redshift z is the time variable
    const auto& z = sf.get_z();

    // H is the Hubble rate
    GiNaC::ex H = FRW::Hub(z);

    // r is the unit line-of-sight vector to Earth
    auto r_sym = sf.make_symbol("r");
    auto r = sf.make_vector(r_sym);
    sf.declare_parameter(r_sym);

    // mu is RSD parameter = r.\hat{k} = r.k / |k|
    auto mu = sf.make_symbol("mu");
    sf.declare_parameter(mu);

    // define halo bias parameters
    auto b1 = sf.make_symbol("b1");
    auto b2 = sf.make_symbol("b2");
    auto b3 = sf.make_symbol("b3");
    auto bG2 = sf.make_symbol("bG2");
    auto bG3 = sf.make_symbol("bG3");
    auto bdG2 = sf.make_symbol("bdG2");
    auto bGamma3 = sf.make_symbol("bGamma3");

    sf.declare_parameter(b1);
    sf.declare_parameter(b2);
    sf.declare_parameter(b3);
    sf.declare_parameter(bG2);
    sf.declare_parameter(bG3);
    sf.declare_parameter(bdG2);
    sf.declare_parameter(bGamma3);

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
    // note that there is no need to symmetrize explicitly; kernels are symmetrized
    // automatically before insertion into the kernel group
    kernel qs_base{iv_qs, sf};
    delta.add(SPT::DA(z) * alpha(q, s, qs_base, sf));
    delta.add(SPT::DB(z) * gamma(q, s, qs_base, sf));

    // third order
    kernel qst_base{iv_qst, sf};
    delta.add((SPT::DD(z) - SPT::DJ(z)) * 2*gamma_bar(s+t, q, alpha_bar(s, t, qst_base, sf), sf));
    delta.add(SPT::DE(z)                * 2*gamma_bar(s+t, q, gamma_bar(s, t, qst_base, sf), sf));
    delta.add((SPT::DF(z) + SPT::DJ(z)) * 2*alpha_bar(s+t, q, alpha_bar(s, t, qst_base, sf), sf));
    delta.add(SPT::DG(z)                * 2*alpha_bar(s+t, q, gamma_bar(s, t, qst_base, sf), sf));
    delta.add(SPT::DJ(z)                * (alpha(s+t, q, gamma_bar(s, t, qst_base, sf), sf)
                                           - 2*alpha(s+t, q, alpha_bar(s, t, qst_base, sf), sf)));

    // compute kernels for the dark matter velocity potential \phi, v = grad phi -> v(k) = i k phi
    auto delta1 = delta.order(1);
    auto delta2 = delta.order(2);
    auto delta3 = delta.order(3);

    auto phi1 = InverseLaplacian(-diff_t(delta1));
    auto phi2 = InverseLaplacian(-diff_t(delta2) - delta1*Laplacian(phi1) - gradgrad(phi1, delta1));
    auto phi3 = InverseLaplacian(-diff_t(delta3)
                                 - delta1*Laplacian(phi2) - delta2*Laplacian(phi1)
                                 - gradgrad(phi1, delta2) - gradgrad(phi2, delta1));

    auto phi = phi1 + phi2 + phi3;


    // build halo overdensity \delta
    // first, need velocity potentials for the Galileon terms
    auto Phi_delta = InverseLaplacian(delta);
    auto Phi_v = -phi/H;

    auto G2 = Galileon2(Phi_delta);
    auto G3 = Galileon3(Phi_delta);
    auto Gamma3 = Galileon2(Phi_delta) - Galileon2(Phi_v);

    auto deltah = b1*delta + (b2/2)*delta*delta + (b3/6)*delta*delta*delta
                  + bG2*G2 + bdG2*G2*delta + bG3*G3 + bGamma3*Gamma3;


    // set up momentum label k, corresponding to external momentum in 2pf
    auto k = sf.make_symbol("k");
    sf.declare_parameter(k);


    // build expression for the redshift-space overdensities,
    // for both dark matter and halos

    // utility function to perform the redshift-space transformation
    // note that we don't have to adjust r_dot_v for the halo power spectrum, so it can just be
    // captured from the exterior scope.
    // Of course, delta has to be adjusted.
    auto r_dot_v = dotgrad(r, phi);
    auto make_delta_rsd = [&](const auto& kmu, const auto& d) -> auto
      {
        return d
               - (GiNaC::I / H) * kmu * r_dot_v
               - (GiNaC::I / H) * kmu * (r_dot_v * d)
               - (GiNaC::numeric{1} / (2*H*H)) * kmu*kmu * (r_dot_v * r_dot_v)
               - (GiNaC::numeric{1} / (2*H*H)) * kmu*kmu * (r_dot_v * r_dot_v * d)
               + (GiNaC::I / (3*2*H*H*H)) * kmu*kmu*kmu * (r_dot_v * r_dot_v * r_dot_v);
      };

    // dark matter in redshift-space
    auto k1mu = k*mu;
    auto k2mu = -k*mu;
    auto delta_rsd_k1 = make_delta_rsd(k1mu, delta);
    auto delta_rsd_k2 = make_delta_rsd(k2mu, delta);

    // halos in redshift-space
    auto deltah_rsd_k1 = make_delta_rsd(k1mu, deltah);
    auto deltah_rsd_k2 = make_delta_rsd(k2mu, deltah);

    // construct 1-loop \delta power spectrum
    Pk_one_loop Pk_delta{deltah_rsd_k1, deltah_rsd_k2, k, sf};

    // simplify mu-dependence
    Pk_delta.canonicalize_external_momenta();
    Pk_delta.simplify(GiNaC::exmap{ {Angular::Cos(k,r_sym), mu} });

    // remove unwanted r factors, which are equal to unity (r is a unit vector)
    Pk_delta.simplify(GiNaC::exmap{ {r_sym, GiNaC::ex{1}} });

//    auto& tree = Pk_delta.get_tree();
//    std::cout << "Tree-level P(k):" << '\n';
//    std::cout << tree << '\n';

//    auto& P13 = Pk_delta.get_13();
//    std::cout << "Loop-level 13 P(k):" << '\n';
//    std::cout << P13 << '\n';

//    auto& P22 = Pk_delta.get_22();
//    std::cout << "Loop-level 22 P(k):" << '\n';
//    std::cout << P22 << '\n';

    // break result into powers of mu, grouped by the bias coefficients involved
    GiNaC_symbol_set filter_syms{b1, b2, b3, bG2, bdG2, bG3, bGamma3};

    Pk_rsd Pk_nobias{Pk_delta, mu, filter_list{}, filter_syms};

    Pk_rsd Pk_b1{Pk_delta, mu, filter_list{ {b1,1} }, filter_syms};
    Pk_rsd Pk_b2{Pk_delta, mu, filter_list{ {b2,1} }, filter_syms};
    Pk_rsd Pk_b3{Pk_delta, mu, filter_list{ {b3,1} }, filter_syms};
    Pk_rsd Pk_bG2{Pk_delta, mu, filter_list{ {bG2,1} }, filter_syms};
    Pk_rsd Pk_bG3{Pk_delta, mu, filter_list{ {bG3,1} }, filter_syms};                           // zero
    Pk_rsd Pk_bdG2{Pk_delta, mu, filter_list{ {bdG2,1} }, filter_syms};
    Pk_rsd Pk_bGamma3{Pk_delta, mu, filter_list{ {bGamma3,1} }, filter_syms};

    Pk_rsd Pk_b1b1{Pk_delta, mu, filter_list{ {b1,2} }, filter_syms};
    Pk_rsd Pk_b1b2{Pk_delta, mu, filter_list{ {b1,1}, {b2,1} }, filter_syms};
    Pk_rsd Pk_b1b3{Pk_delta, mu, filter_list{ {b1,1}, {b3,1} }, filter_syms};
    Pk_rsd Pk_b2b2{Pk_delta, mu, filter_list{ {b2,2} }, filter_syms};

    Pk_rsd Pk_b1bG2{Pk_delta, mu, filter_list{ {b1,1}, {bG2,1} }, filter_syms};
    Pk_rsd Pk_bG2bG2{Pk_delta, mu, filter_list{ {bG2,2} }, filter_syms};
    Pk_rsd Pk_b2bG2{Pk_delta, mu, filter_list{ {b2,1}, {bG2,1} }, filter_syms};
    Pk_rsd Pk_b1bG3{Pk_delta, mu, filter_list{ {b1,1}, {bG3,1} }, filter_syms};                 // zero
    Pk_rsd Pk_b1bdG2{Pk_delta, mu, filter_list{ {b1,1}, {bdG2,1} }, filter_syms};
    Pk_rsd Pk_b1bGamma3{Pk_delta, mu, filter_list{ {b1,1}, {bGamma3,1} }, filter_syms};

    std::cout << "** COUNTERTERM MAP" << '\n' << '\n';

    std::cout << "OPERATOR MIXING:" << '\n';
    write_map(Pk_nobias, Pk_b1, Pk_b2, Pk_b3, Pk_bG2, Pk_bdG2, Pk_bGamma3,
              Pk_b1b1, Pk_b1b2, Pk_b1b3, Pk_b2b2,
              Pk_b1bG2, Pk_bG2bG2, Pk_b2bG2, Pk_b1bdG2,
              Pk_b1bGamma3, k, mixing_divergences);

    std::cout << "STOCHASTIC COUNTERTERMS:" << '\n';
    write_map(Pk_nobias, Pk_b1, Pk_b2, Pk_b3, Pk_bG2, Pk_bdG2, Pk_bGamma3,
              Pk_b1b1, Pk_b1b2, Pk_b1b3, Pk_b2b2,
              Pk_b1bG2, Pk_bG2bG2, Pk_b2bG2, Pk_b1bdG2,
              Pk_b1bGamma3, k, stochastic_divergences);

//    std::cout << '\n' << '\n';
//
//    std::cout << "** NUMBER OF KERNELS" << '\n' << '\n';
//    count_kernels(Pk_nobias, Pk_b1, Pk_b2, Pk_b3, Pk_bG2, Pk_bdG2, Pk_bGamma3,
//                  Pk_b1b1, Pk_b1b2, Pk_b1b3, Pk_b2b2,
//                  Pk_b1bG2, Pk_bG2bG2, Pk_b2bG2, Pk_b1bdG2,
//                  Pk_b1bGamma3);

    return EXIT_SUCCESS;
  }
