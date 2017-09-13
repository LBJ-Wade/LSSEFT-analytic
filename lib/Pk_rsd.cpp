//
// Created by David Seery on 11/09/2017.
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

#include "Pk_rsd.h"


Pk_rsd_group::Pk_rsd_group(GiNaC::symbol mu_, filter_list pt_, GiNaC_symbol_set sy_)
  : mu(std::move(mu_)),
    pattern(std::move(pt_)),
    filter_symbols(std::move(sy_))
  {
    // set up an exmap that will set all filter_symbols to zero; we use this after pattern matching
    // to kill any remaining terms; eg. asking for the term linear in b1 from b1(1+b2) will give 1+b2, and we
    // need to set b2 to zero to get the required result '1'
    for(const auto& sym : filter_symbols)
      {
        filter_map[sym] = 0;
      }
  }


void Pk_rsd_group::emplace(const one_loop_element& elt)
  {
    auto filter = [&](const auto& e, unsigned int order) -> auto
      {
        // copy input element; one_loop_element doesn't capture any of its data by reference,
        // so what we get is a standalone copy.
        // That means we can later modify the Pk_one_loop object if we require
        auto f = std::make_unique<one_loop_element>(e);

        // filter for pattern
        for(const auto& pat : this->pattern)
          {
            f->filter(pat.first, pat.second);
          }

        // kill any remaining symbols (eg. if we are looking for the term linear in b1 we don't get b1 b2 contributions)
        f->simplify(this->filter_map);

        // filter for power of mu
        f->filter(this->mu, order);

        return f;
      };

    auto mu0 = filter(elt, 0);
    auto mu2 = filter(elt, 2);
    auto mu4 = filter(elt, 4);
    auto mu6 = filter(elt, 6);
    auto mu8 = filter(elt, 8);

    one_loop_element_key mu0_key{*mu0};
    one_loop_element_key mu2_key{*mu2};
    one_loop_element_key mu4_key{*mu4};
    one_loop_element_key mu6_key{*mu6};
    one_loop_element_key mu8_key{*mu8};

    if(!mu0->null()) this->mu0[mu0_key].push_back(std::move(mu0));
    if(!mu2->null()) this->mu2[mu2_key].push_back(std::move(mu2));
    if(!mu4->null()) this->mu4[mu4_key].push_back(std::move(mu4));
    if(!mu6->null()) this->mu6[mu6_key].push_back(std::move(mu6));
    if(!mu8->null()) this->mu8[mu8_key].push_back(std::move(mu8));
  }


GiNaC::exvector Pk_rsd_group::get_UV_limit(unsigned int order) const
  {
    GiNaC::exvector values;

    auto build = [&](const auto& db) -> auto
      {
        GiNaC::ex value{0};

        for(const auto& item : db)
          {
            for(const auto& iptr : item.second)
              {
                const auto& i = *iptr;
                value += i.get_UV_limit(order);
              }
          }

        return value;
      };

    values.push_back(build(this->mu0));
    values.push_back(build(this->mu2));
    values.push_back(build(this->mu4));
    values.push_back(build(this->mu6));
    values.push_back(build(this->mu8));

    return values;
  }


void Pk_rsd_group::write(std::ostream& out) const
  {
    out << "-- mu^0" << '\n'; if(!this->mu0.empty()) out << this->mu0 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^2" << '\n'; if(!this->mu0.empty()) out << this->mu2 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^4" << '\n'; if(!this->mu0.empty()) out << this->mu4 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^6" << '\n'; if(!this->mu0.empty()) out << this->mu6 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^8" << '\n'; if(!this->mu0.empty()) out << this->mu8 << '\n'; else out << "   <empty>" << '\n';
  }


std::vector<size_t> Pk_rsd_group::get_number_time_functions() const
  {
    std::vector<size_t> values;

    values.push_back(this->mu0.size());
    values.push_back(this->mu2.size());
    values.push_back(this->mu4.size());
    values.push_back(this->mu6.size());
    values.push_back(this->mu8.size());

    return values;
  }


Pk_rsd::Pk_rsd(const Pk_one_loop& Pk, const GiNaC::symbol& mu_, const filter_list pt_, const GiNaC_symbol_set sy_)
  : mu(mu_),
    pattern(pt_),
    Ptree(mu_, pt_, sy_),
    P13(mu_, pt_, sy_),
    P22(mu_, pt_, sy_)
  {
    // extract Pk_db elements from Pk
    this->filter(Ptree, Pk.get_tree());
    this->filter(P13, Pk.get_13());
    this->filter(P22, Pk.get_22());
  }


void Pk_rsd::filter(Pk_rsd_group& dest, const Pk_one_loop_impl::Pk_db& source)
  {
    for(const auto& item : source)
      {
        if(!item.second) continue;    // skip if pointer is empty

        const auto& reduced = *item.second;
        const auto& db  = reduced.get_db();

        for(const auto& term : db)
          {
            for(const auto& t : term.second)
              {
                const auto& elt = *t;
                dest.emplace(elt);
              }
          }
      }
  }


void Pk_rsd::write(std::ostream& out) const
  {
    out << "Tree-level:" << '\n';
    out << this->Ptree << '\n';
    out << "Loop-level 13 terms:" << '\n';
    out << this->P13 << '\n';
    out << "Loop-level 22 terms:" << '\n';
    out << this->P22 << '\n';
  }


std::ostream& operator<<(std::ostream& str, const Pk_rsd& obj)
  {
    obj.write(str);
    return str;
  }


std::ostream& operator<<(std::ostream& str, const Pk_rsd_group& obj)
  {
    obj.write(str);
    return str;
  }
