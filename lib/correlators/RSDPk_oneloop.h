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

#ifndef LSSEFT_ANALYTIC_PK_RSD_H
#define LSSEFT_ANALYTIC_PK_RSD_H


#include "RSDPk_set.h"


class RSDPk_oneloop
  {

    // TYPES

  public:

    //! tree-level RSD Pk set
    using RSDPk_tree_set = RSDPk_set<oneloop_element_db>;

    //! one-loop RSD Pk set
    using RSDPk_oneloop_set = RSDPk_set<oneloop_element_db>;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor accepts a Pk_one_loop and decomposes it into powers of mu.
    //! filters terms according to the provided pattern
    RSDPk_oneloop(const Pk_oneloop& Pk, const GiNaC::symbol& mu_, const filter_list pt_, const GiNaC_symbol_set sy_,
                  bool v=false);

    //! destructor is default
    ~RSDPk_oneloop() = default;


    // INTERNAL API

  protected:

    //! filter a Pk_db into a destination Pk_rsd_group
    template <typename SetType>
    void filter(SetType& dest, const Pk_oneloop::Pk_db& source);


    // ACCESSORS

  public:

    //! get tree-level terms
    const RSDPk_tree_set& get_tree() const { return this->Ptree; }

    //! get 13 terms
    const RSDPk_oneloop_set& get_13() const { return this->P13; }

    //! get 22 terms
    const RSDPk_oneloop_set& get_22() const { return this->P22; }


    // SERVICES

  public:

    //! write self to stream
    void write(std::ostream& out) const;

    //! write Mathematica script for 13 and 22 integrals at each power of mu
    void write_Mathematica(std::ostream& out) const;


    // INTERNAL DATA

  private:

    //! cache mu parameter
    const GiNaC::symbol mu;

    //! cache tag used to represent this power spectrum group
    std::string tag;


    // FILTER/PATTERN

    //! cache pattern used to filter input expression
    const filter_list pattern;

    //! cache GiNaC expression giving symbolic filter list
    GiNaC::ex symbolic_filter;


    // DATABASES

    //! tree contributions
    RSDPk_tree_set Ptree;

    //! 13 contributions
    RSDPk_oneloop_set P13;

    //! 22 contributions
    RSDPk_oneloop_set P22;

  };


template <typename SetType>
void RSDPk_oneloop::filter(SetType& dest, const Pk_oneloop::Pk_db& source)
  {
    // walk through the source database, filtering contributions to the reduced integral (if present)
    // that match our configured symbol set
    // Matching contributions get pushed into the destination database 'dest',
    // which will be either tree, 13 or 22

    for(const auto& item : source)
      {
        const Pk_oneloop::Pk_db::expr_type& lp = *item.second.first;
        const Pk_oneloop::Pk_db::reduced_ptr_type& ri = item.second.second;

        if(!ri) continue;    // skip if pointer is empty

        // get database of one-loop-reduced-integral elements
        const auto& db  = ri->get_db();

        // walk through this list
        for(const auto& record : db)
          {
            const std::unique_ptr<oneloop_element>& elt = record.second;
            if(elt) dest.emplace(*elt);
          }
      }
  }


//! stream insertion operator
std::ostream& operator<<(std::ostream& str, const RSDPk_oneloop& obj);


#endif //LSSEFT_ANALYTIC_PK_RSD_H
