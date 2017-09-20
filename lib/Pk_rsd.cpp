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

    this->emplace(filter(elt, 0), this->mu0);
    this->emplace(filter(elt, 2), this->mu2);
    this->emplace(filter(elt, 4), this->mu4);
    this->emplace(filter(elt, 6), this->mu6);
    this->emplace(filter(elt, 8), this->mu8);
  }


GiNaC::exvector Pk_rsd_group::get_UV_limit(unsigned int order) const
  {
    GiNaC::exvector values;

    auto build = [&](const one_loop_element_db& db) -> auto
      {
        GiNaC::ex value{0};

        for(const auto& record : db)
          {
            const auto& data = record.second;

            if(data) value += data->get_UV_limit(order);
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
    out << "-- mu^2" << '\n'; if(!this->mu2.empty()) out << this->mu2 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^4" << '\n'; if(!this->mu4.empty()) out << this->mu4 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^6" << '\n'; if(!this->mu6.empty()) out << this->mu6 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^8" << '\n'; if(!this->mu8.empty()) out << this->mu8 << '\n'; else out << "   <empty>" << '\n';
  }


std::vector< std::vector<time_function> > Pk_rsd_group::get_time_functions() const
  {
    std::vector< std::vector<time_function> > values(5);

    auto build = [](const one_loop_element_db& db, std::vector<time_function>& dest)
      {
        for(const auto& item : db)
          {
            const GiNaC::ex& tm = item.first.get_time_function();

            auto it = std::find_if(dest.begin(), dest.end(),
                                   [&](const time_function& t) -> bool
                                     { return static_cast<bool>(t == tm); });

            if(it == dest.end()) dest.push_back(tm);
          }
      };

    build(this->mu0, values[0]);
    build(this->mu2, values[1]);
    build(this->mu4, values[2]);
    build(this->mu6, values[3]);
    build(this->mu8, values[4]);

    return values;
  }


void Pk_rsd_group::emplace(std::unique_ptr<one_loop_element> elt, one_loop_element_db& db)
  {
    one_loop_element_key key{*elt};

    // check whether a compatible entry is already present in the database
    auto it = db.find(key);


    // if such an entry is present we can just compose the integrands
    if(it != db.end())
      {
        *it->second += *elt;
        return;
      }

    // otherwise, we need to insert a new element
    auto res = db.emplace(std::move(key), std::move(elt));
    if(!res.second) throw exception(ERROR_ONE_LOOP_ELEMENT_INSERT_FAILED_RSD, exception_code::loop_integral_error);
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
        const loop_integral& lp = *item.second.first;
        const std::unique_ptr<one_loop_reduced_integral>& ri = item.second.second;

        if(!ri) continue;    // skip if pointer is empty

        const auto& db  = ri->get_db();

        for(const auto& record : db)
          {
            const std::unique_ptr<one_loop_element>& elt = record.second;
            if(elt) dest.emplace(*elt);
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
