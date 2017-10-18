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
#include <fstream>
#include <functional>

#include "services/service_locator.h"

#include "lib/vector.h"
#include "lib/initial_value.h"
#include "lib/fourier_kernel.h"
#include "lib/Pk_one_loop.h"
#include "lib/Pk_rsd.h"
#include "lib/detail/special_functions.h"

#include "SPT/time_functions.h"
#include "SPT/one_loop_kernels.h"

#include "backends/LSSEFT.h"

#include "instruments/timing_instrument.h"


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
    return generate_map(rsd_22, k, 4, 1);
  }


std::vector<std::string> mixing_divergences(const Pk_rsd& rsd, const GiNaC::symbol& k)
  {
    const auto& rsd_13 = rsd.get_13();
    return generate_UV_limit(rsd_13, k, 4, 1);
  }


std::vector<std::string> stochastic_divergences(const Pk_rsd& rsd, const GiNaC::symbol& k)
  {
    const auto& rsd_22 = rsd.get_22();
    return generate_UV_limit(rsd_22, k, 4, 1);
  }


template <typename MapGenerator>
void write_map(const Pk_rsd_set& Pks, const GiNaC::symbol& k, MapGenerator make_map)
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

    for(const auto& item : Pks)
      {
        const std::string& name = item.first;
        const Pk_rsd& Pk = item.second.get();

        write(name, Pk);
      }
  }


int main(int argc, char* argv[])
  {
    // generate service objects
    symbol_factory sf;
    argument_cache args{argc, argv};

    // build service locator
    service_locator loc{args, sf};

    if(!args.get_counterterms() && args.get_output_path().empty() && args.get_Mathematica_output().empty())
      exit(EXIT_SUCCESS);

    // redshift z is the time variable
    const auto& z = sf.get_z();

    // H is the Hubble rate, f is linear growth factor
    GiNaC::ex H = FRW::Hub(z);
    GiNaC::ex f = SPT::f(z);

    // r is the unit line-of-sight vector to Earth
    auto r_sym = sf.make_symbol("r");
    auto r = sf.make_vector(r_sym);
    sf.declare_parameter(r_sym);

    // mu is RSD parameter = r.\hat{k} = r.k / |k|
    auto mu = sf.make_symbol("mu");
    sf.declare_parameter(mu);

    // define halo bias parameters
    auto b1_1 = sf.make_symbol("b1_1");
    auto b1_2 = sf.make_symbol("b1_2");
    auto b1_3 = sf.make_symbol("b1_3");

    auto b2_2 = sf.make_symbol("b2_2");
    auto b2_3 = sf.make_symbol("b2_3");

    auto bG2_2 = sf.make_symbol("bG2_2");
    auto bG2_3 = sf.make_symbol("bG2_3");

    auto b3 = sf.make_symbol("b3");
    auto bG3 = sf.make_symbol("bG3");
    auto bdG2 = sf.make_symbol("bdG2");
    auto bGamma3 = sf.make_symbol("bGamma3");

    sf.declare_parameter(b1_1).declare_parameter(b1_2).declare_parameter(b1_3);
    sf.declare_parameter(b2_2).declare_parameter(b2_3);
    sf.declare_parameter(bG2_2).declare_parameter(bG2_3);
    sf.declare_parameter(b3);
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


    auto delta_timer = std::make_unique<timing_instrument>("Construct \\delta Fourier representation");

    // set up kernels for the dark matter overdensity \delta
    auto delta = loc.make_fourier_kernel<3>();

    // linear order
    delta.add(SPT::D(z), iv_q, 1);

    // second order
    // we don't symmetrize explicitly; kernels are symmetrized automatically
    // if this feature is not disabled
    kernel qs_base{iv_qs, loc};
    delta.add(SPT::DA(z) * alpha(q, s, qs_base, loc));
    delta.add(SPT::DB(z) * gamma(q, s, qs_base, loc));

    // third order
    kernel qst_base{iv_qst, loc};
    delta.add((SPT::DD(z) - SPT::DJ(z)) * 2*gamma_bar(s+t, q, alpha_bar(s, t, qst_base, loc), loc));
    delta.add(SPT::DE(z)                * 2*gamma_bar(s+t, q, gamma_bar(s, t, qst_base, loc), loc));
    delta.add((SPT::DF(z) + SPT::DJ(z)) * 2*alpha_bar(s+t, q, alpha_bar(s, t, qst_base, loc), loc));
    delta.add(SPT::DG(z)                * 2*alpha_bar(s+t, q, gamma_bar(s, t, qst_base, loc), loc));
    delta.add(SPT::DJ(z)                * (alpha(s+t, q, gamma_bar(s, t, qst_base, loc), loc)
                                           - 2*alpha(s+t, q, alpha_bar(s, t, qst_base, loc), loc)));

    // extract different orders of \delta
    auto delta_1 = delta.order(1);
    auto delta_2 = delta.order(2);
    auto delta_3 = delta.order(3);

    // extract different orders of \delta^2
    auto deltasq = delta*delta;
    auto deltasq_2 = deltasq.order(2);
    auto deltasq_3 = deltasq.order(3);

    delta_timer.reset(nullptr);


    auto phi_timer = std::make_unique<timing_instrument>("Construct velocity potential \\phi");

    // compute kernels for the dark matter velocity potential \phi, v = grad phi -> v(k) = i k phi
    auto phi1 = InverseLaplacian(-diff_t(delta_1));
    auto phi2 = InverseLaplacian(-diff_t(delta_2) - delta_1*Laplacian(phi1) - gradgrad(phi1, delta_1));
    auto phi3 = InverseLaplacian(-diff_t(delta_3)
                                 - delta_1*Laplacian(phi2) - delta_2*Laplacian(phi1)
                                 - gradgrad(phi1, delta_2) - gradgrad(phi2, delta_1));

    auto phi = phi1 + phi2 + phi3;

    phi_timer.reset(nullptr);


    auto Galileon_timer = std::make_unique<timing_instrument>("Construct Galileon operators");

    // build halo overdensity \delta
    // first, need velocity potentials for the Galileon terms
    auto Phi_delta = InverseLaplacian(delta);
    auto Phi_v = -phi/H;

    auto G2 = Galileon2(Phi_delta);
    auto G2_2 = G2.order(2);
    auto G2_3 = G2.order(3);

    auto G3 = Galileon3(Phi_delta);
    auto Gamma3 = Galileon2(Phi_delta) - Galileon2(Phi_v);

    Galileon_timer.reset(nullptr);


    auto halo_timer = std::make_unique<timing_instrument>("Construct halo overdensity field");

    auto vp1 = phi1 / (H*f);
    auto vp2 = phi2 / (H*f);

    auto deltah_b1 = b1_1*delta_1 + b1_2*delta_2 + b1_3*delta_3;

    auto deltah_b1_adv = - (b1_1 - b1_2) * gradgrad(vp1, delta_1)
                         - (b1_2 - b1_3) * gradgrad(vp1, delta_2)
                         - (b1_1 - b1_3) * gradgrad(vp2, delta_1) / 2
                         + ((b1_1 + b1_3) / 2 - b1_2) * convective_bias_term(vp1, delta_1);

    auto deltah_b2 = (b2_2/2)*deltasq_2 + (b2_3/2)*deltasq_3;

    auto deltah_b2_adv = - (b2_2/2-b2_3/2) * gradgrad(vp1, deltasq_2);

    auto deltah_G2 = bG2_2*G2_2 + bG2_3*G2_3;

    auto deltah_G2_adv = - (bG2_2 - bG2_3) * gradgrad(vp1, G2_2);

    auto deltah_cubic = (b3/6)*delta*delta*delta + bdG2*G2*delta + bG3*G3 + bGamma3*Gamma3;

    auto deltah = deltah_b1 + deltah_b1_adv + deltah_b2 + deltah_b2_adv + deltah_G2 + deltah_G2_adv + deltah_cubic;

    halo_timer.reset(nullptr);


    // set up momentum label k, corresponding to external momentum in 2pf
    auto k = sf.make_symbol("k");
    sf.declare_parameter(k);


    Pk_one_loop test_Pk_b1{"delta_1 x delta_3", "d", delta_1,
                           delta_3 + gradgrad(vp1, delta_2) + gradgrad(vp2, delta_1) / 2 +
                           convective_bias_term(vp1, delta_1) / 2, k, loc};
    Pk_one_loop test_Pk_b2{"delta_1 x delta^2_3", "dd", delta_1, deltasq_3 + gradgrad(vp1, deltasq_2), k, loc};
    Pk_one_loop test_Pk_b3{"delta_1 x delta^3", "ddd", delta_1, delta*delta*delta, k, loc};
    Pk_one_loop test_Pk_bG2{"delta_1 x G2_3", "G2", delta_1, G2_3 + gradgrad(vp1, G2_2), k, loc};
    Pk_one_loop test_Pk_bdG2{"delta_1 x G2 delta", "G2d", delta_1, G2*delta, k, loc};
    Pk_one_loop test_Pk_Gamma3{"delta_1 x Gamma3", "Gamma3", delta_1, Gamma3, k, loc};

    if(!args.get_Mathematica_output().empty())
      {
        std::ofstream mma_out{args.get_Mathematica_output().string(), std::ios_base::out | std::ios_base::trunc};
        test_Pk_b1.write_Mathematica(mma_out);
        test_Pk_b2.write_Mathematica(mma_out);
        test_Pk_b3.write_Mathematica(mma_out);
        test_Pk_bG2.write_Mathematica(mma_out);
        test_Pk_bdG2.write_Mathematica(mma_out);
        test_Pk_Gamma3.write_Mathematica(mma_out);
        mma_out.close();
      }


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

    auto deltarsd_timer = std::make_unique<timing_instrument>("RSD transform for dark matter overdensity");
    auto delta_rsd_k1 = make_delta_rsd(k1mu, delta);
    auto delta_rsd_k2 = make_delta_rsd(k2mu, delta);
    deltarsd_timer.reset(nullptr);

    // halos in redshift-space
    auto deltahrsd_timer = std::make_unique<timing_instrument>("RSD transform for halo overdensity");
    auto deltah_rsd_k1 = make_delta_rsd(k1mu, deltah);
    auto deltah_rsd_k2 = make_delta_rsd(k2mu, deltah);
    deltahrsd_timer.reset(nullptr);

    // construct 1-loop \delta power spectrum
    auto Pk_timer = std::make_unique<timing_instrument>("Construct 1-loop power spectrum");
    Pk_one_loop Pk_delta{"1-loop halo RSD P(k)", "halo", deltah_rsd_k1, deltah_rsd_k2, k, loc};

    // simplify mu-dependence
    Pk_delta.canonicalize_external_momenta();
    Pk_delta.simplify(GiNaC::exmap{ {Angular::Cos(k,r_sym), mu} });

    // remove unwanted r factors, which are equal to unity (r is a unit vector)
    Pk_delta.simplify(GiNaC::exmap{ {r_sym, GiNaC::ex{1}} });

//    if(!args.get_Mathematica_output().empty())
//      {
//        std::ofstream mma_out{args.get_Mathematica_output().string(), std::ios_base::out | std::ios_base::trunc};
//        Pk_delta.write_Mathematica(mma_out);
//        mma_out.close();
//      }

    Pk_timer.reset(nullptr);

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
    GiNaC_symbol_set filter_syms{b1_1, b1_2, b1_3, b2_2, b2_3, b3, bG2_2, bG2_3, bdG2, bG3, bGamma3};

    auto Pk_rsd_timer = std::make_unique<timing_instrument>("Extract RSD mu coefficients");

    Pk_rsd Pk_nobias{Pk_delta, mu, filter_list{}, filter_syms};

    Pk_rsd Pk_b1_1{Pk_delta, mu, filter_list{ {b1_1,1} }, filter_syms};
    Pk_rsd Pk_b1_2{Pk_delta, mu, filter_list{ {b1_2,1} }, filter_syms};
    Pk_rsd Pk_b1_3{Pk_delta, mu, filter_list{ {b1_3,1} }, filter_syms};

    Pk_rsd Pk_b2_2{Pk_delta, mu, filter_list{ {b2_2,1} }, filter_syms};
    Pk_rsd Pk_b2_3{Pk_delta, mu, filter_list{ {b2_3,1} }, filter_syms};

    Pk_rsd Pk_bG2_2{Pk_delta, mu, filter_list{ {bG2_2,1} }, filter_syms};
    Pk_rsd Pk_bG2_3{Pk_delta, mu, filter_list{ {bG2_3,1} }, filter_syms};

    Pk_rsd Pk_b3{Pk_delta, mu, filter_list{ {b3,1} }, filter_syms};
    Pk_rsd Pk_bG3{Pk_delta, mu, filter_list{ {bG3,1} }, filter_syms};                           // zero
    Pk_rsd Pk_bdG2{Pk_delta, mu, filter_list{ {bdG2,1} }, filter_syms};
    Pk_rsd Pk_bGamma3{Pk_delta, mu, filter_list{ {bGamma3,1} }, filter_syms};

    Pk_rsd Pk_b1_1_b1_1{Pk_delta, mu, filter_list{ {b1_1,2} }, filter_syms};
    Pk_rsd Pk_b1_2_b1_2{Pk_delta, mu, filter_list{ {b1_2,2} }, filter_syms};
    Pk_rsd Pk_b1_1_b1_2{Pk_delta, mu, filter_list{ {b1_1,1}, {b1_2,1} }, filter_syms};
    Pk_rsd Pk_b1_1_b1_3{Pk_delta, mu, filter_list{ {b1_1,1}, {b1_3,1} }, filter_syms};

    Pk_rsd Pk_b1_1_b2_2{Pk_delta, mu, filter_list{ {b1_1,1}, {b2_2,1} }, filter_syms};
    Pk_rsd Pk_b1_1_b2_3{Pk_delta, mu, filter_list{ {b1_1,1}, {b2_3,1} }, filter_syms};
    Pk_rsd Pk_b1_2_b2_2{Pk_delta, mu, filter_list{ {b1_2,1}, {b2_2,1} }, filter_syms};

    Pk_rsd Pk_b1_1_b3{Pk_delta, mu, filter_list{ {b1_1,1}, {b3,1} }, filter_syms};

    Pk_rsd Pk_b2_2_b2_2{Pk_delta, mu, filter_list{ {b2_2,2} }, filter_syms};

    Pk_rsd Pk_b1_1_bG2_2{Pk_delta, mu, filter_list{ {b1_1,1}, {bG2_2,1} }, filter_syms};
    Pk_rsd Pk_b1_1_bG2_3{Pk_delta, mu, filter_list{ {b1_1,1}, {bG2_3,1} }, filter_syms};
    Pk_rsd Pk_b1_2_bG2_2{Pk_delta, mu, filter_list{ {b1_2,1}, {bG2_2,1} }, filter_syms};

    Pk_rsd Pk_bG2_2_bG2_2{Pk_delta, mu, filter_list{ {bG2_2,2} }, filter_syms};

    Pk_rsd Pk_b2_2_bG2_2{Pk_delta, mu, filter_list{ {b2_2,1}, {bG2_2,1} }, filter_syms};

    Pk_rsd Pk_b1_1_bG3{Pk_delta, mu, filter_list{ {b1_1,1}, {bG3,1} }, filter_syms};            // zero

    Pk_rsd Pk_b1_1_bdG2{Pk_delta, mu, filter_list{ {b1_1,1}, {bdG2,1} }, filter_syms};

    Pk_rsd Pk_b1_1_bGamma3{Pk_delta, mu, filter_list{ {b1_1,1}, {bGamma3,1} }, filter_syms};

    Pk_rsd_timer.reset(nullptr);


    Pk_rsd_set Pks =
      {
        {"nobias", std::ref(Pk_nobias)},
        {"b1_1", std::ref(Pk_b1_1)},
        {"b1_2", std::ref(Pk_b1_2)},
        {"b1_3", std::ref(Pk_b1_3)},
        {"b2_2", std::ref(Pk_b2_2)},
        {"b2_3", std::ref(Pk_b2_3)},
        {"bG2_2", std::ref(Pk_bG2_2)},
        {"bG2_3", std::ref(Pk_bG2_3)},
        {"b3", std::ref(Pk_b3)},
        {"bdG2", std::ref(Pk_bdG2)},
        {"bGamma3", std::ref(Pk_bGamma3)},
        {"b1_1_b1_1", std::ref(Pk_b1_1_b1_1)},
        {"b1_2_b1_2", std::ref(Pk_b1_2_b1_2)},
        {"b1_1_b1_2", std::ref(Pk_b1_1_b1_2)},
        {"b1_1_b1_3", std::ref(Pk_b1_1_b1_3)},
        {"b1_1_b2_2", std::ref(Pk_b1_1_b2_2)},
        {"b1_1_b2_3", std::ref(Pk_b1_1_b2_3)},
        {"b1_2_b2_2", std::ref(Pk_b1_2_b2_2)},
        {"b1_1_b3", std::ref(Pk_b1_1_b3)},
        {"b1_1_bG2_2", std::ref(Pk_b1_1_bG2_2)},
        {"b1_1_bG2_3", std::ref(Pk_b1_1_bG2_3)},
        {"b1_2_bG2_2", std::ref(Pk_b1_2_bG2_2)},
        {"bG2_2_bG2_2", std::ref(Pk_bG2_2_bG2_2)},
        {"b2_2_bG2_2", std::ref(Pk_b2_2_bG2_2)},
        {"b1_1_bdG2", std::ref(Pk_b1_1_bdG2)},
        {"b1_1_bGamma3", std::ref(Pk_b1_1_bGamma3)},
      };

    if(args.get_counterterms())
      {
        std::cout << "** COUNTERTERM MAP" << '\n' << '\n';


        std::cout << "OPERATOR MIXING:" << '\n';
        write_map(Pks, k, mixing_divergences);

        std::cout << "STOCHASTIC COUNTERTERMS:" << '\n';
        write_map(Pks, k, stochastic_divergences);
      }

    if(!args.get_output_path().empty())
      {
        LSSEFT backend{args.get_output_path(), loc};
        backend.add(Pks);

        backend.write();
      }


    return EXIT_SUCCESS;
  }
