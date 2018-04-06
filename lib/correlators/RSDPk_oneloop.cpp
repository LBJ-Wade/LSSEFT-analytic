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

#include <algorithm>

#include "RSDPk_oneloop.h"


RSDPk_oneloop::RSDPk_oneloop(const Pk_oneloop& Pk, const GiNaC::symbol& mu_,
                             const filter_list pt_, const GiNaC_symbol_set sy_, bool v)
  : mu(mu_),
    pattern(pt_),
    Ptree(mu_, pt_, sy_, std::string{"tree"}, v),
    P13(mu_, pt_, sy_, std::string{"13"}, v),
    P22(mu_, pt_, sy_, std::string{"22"}, v),
    symbolic_filter(1)
  {
    // build tag from filter list
    std::ostringstream tag_str;
    tag_str << Pk.get_tag();

    for(const auto& sym : pt_)
      {
        // remove occurrences of '_' in symbol name, since these aren't legal in Mathematica symbols
        // and we use the tag for that purpose
        std::string name = sym.first.get_name();
        std::replace(name.begin(), name.end(), '_', 'z');

        for(unsigned int i = 0; i < sym.second; ++i)
          {
            tag_str << "z" << name;
          }
      }

    tag = tag_str.str();

    // build up symbolic representation of filter pattern
    for(const auto& item : pattern)
      {
        symbolic_filter *= GiNaC::pow(item.first, item.second);
      }

    // filter elements (term by term) from parent Pk_oneloop according to whether they match
    // the specific symbol set
    this->filter(Ptree, Pk.get_tree());
    this->filter(P13, Pk.get_13());
    this->filter(P22, Pk.get_22());

    // remove empty records
    this->Ptree.prune();
    this->P13.prune();
    this->P22.prune();

    error_handler err;
    if(this->Ptree.empty() && this->P13.empty() && this->P22.empty())
      {
        std::ostringstream msg;
        msg << WARNING_PK_RSD_EMPTY << " '" << this->symbolic_filter << "'";
        err.warn(msg.str());
      }
  }


void RSDPk_oneloop::filter(RSDPk_oneloop::RSDPk_tree_set& dest, const Pk_oneloop::Pk_tree_db& source)
  {
    // Walk through the source database, pushing its oneloop_element entries into the destination set.
    // If the destination set has a filter in place then this will automatically apply the filter
    // against any entries we push, so there is no need to check that here

    for(const auto& item : source)
      {
        // extract oneloop_expression object and oneloop_reduced_integral object
        const Pk_oneloop::Pk_tree_db::expr_type& elt = *item.second;

        // push into destination set
        // Filtering occurs automatically if required
        dest.emplace(elt);
      }
  }

void RSDPk_oneloop::filter(RSDPk_oneloop_set& dest, const Pk_oneloop::Pk_oneloop_db& source)
  {
    // Walk through the source database, pushing its oneloop_element entries into the destination set.
    // If the destination set has a filter in place then this will automatically apply the filter
    // against any entries we push, so there is no need to check that here

    for(const auto& item : source)
      {
        // extract oneloop_expression object and oneloop_reduced_integral object
        const Pk_oneloop::Pk_oneloop_db::expr_type& lp = *item.second.first;
        const Pk_oneloop::Pk_oneloop_db::reduced_ptr_type& ri = item.second.second;

        // if no reduced integral, then nothing to insert
        if(!ri) continue;    // skip if pointer is empty

        // get database of one-loop-reduced-integral elements
        const auto& db  = ri->get_db();

        // walk through this list
        for(const auto& record : db)
          {
            const std::unique_ptr<oneloop_element>& elt = record.second;

            // push into destination set
            // Filtering occurs automatically if required
            if(elt) dest.emplace(*elt);
          }
      }
  }


void RSDPk_oneloop::write(std::ostream& out) const
  {
    out << "Tree-level:" << '\n';
    out << this->Ptree << '\n';
    out << "Loop-level 13 terms:" << '\n';
    out << this->P13 << '\n';
    out << "Loop-level 22 terms:" << '\n';
    out << this->P22 << '\n';
  }


void RSDPk_oneloop::write_Mathematica(std::ostream& out) const
  {
    this->Ptree.write_Mathematica(out, this->tag + "Tree", false);
    this->P13.write_Mathematica(out, this->tag + "P13", true);
    this->P22.write_Mathematica(out, this->tag + "P22", false);
  }


std::ostream& operator<<(std::ostream& str, const RSDPk_oneloop& obj)
  {
    obj.write(str);
    return str;
  }
