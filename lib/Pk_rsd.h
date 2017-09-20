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


#include <iostream>

#include "Pk_one_loop.h"


using filter_list = std::vector< std::pair< GiNaC::symbol, unsigned int > >;


class Pk_rsd_group
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    Pk_rsd_group(GiNaC::symbol mu_, filter_list pt_, GiNaC_symbol_set sy_);

    //! destructor
    ~Pk_rsd_group() = default;


    // EMPLACE TERMS

  public:

    //! accepts a one_loop_element and breaks it into individual powers of mu
    void emplace(const one_loop_element& elt);

  protected:

    //! perform emplace on a specific database
    void emplace(std::unique_ptr<one_loop_element> elt, one_loop_element_db& db);


    // TOOLS

  public:

    //! construct UV limit
    GiNaC::exvector get_UV_limit(unsigned int order=2) const;

    //! query number of distinct time functions at each mu
    std::vector< std::vector<time_function> > get_time_functions() const;

    //! prune empty records from the database
    void prune();

  protected:

    //! prune a specific database
    void prune(one_loop_element_db& db);


    // SERVICES

  public:

    //! write self to stream
    void write(std::ostream& out) const;


    // INTERNAL DATA

  private:

    //! cache reference to angular variable mu
    const GiNaC::symbol mu;

    //! cache filtering pattern
    const filter_list pattern;

    //! set of filtering symbols (used to kill any leftover terms)
    const GiNaC_symbol_set filter_symbols;

    //! exmap build from these symbols
    GiNaC::exmap filter_map;


    // TERM-BY-TERM DATABASES

    //! mu^0 terms
    one_loop_element_db mu0;

    //! mu^2 terms
    one_loop_element_db mu2;

    //! mu^4 terms
    one_loop_element_db mu4;

    //! mu^6 terms
    one_loop_element_db mu6;

    //! mu^8 terms;
    one_loop_element_db mu8;

  };


//! stream insertion
std::ostream& operator<<(std::ostream& str, const Pk_rsd_group& obj);


class Pk_rsd
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor accepts a Pk_one_loop and decomposes it into powers of mu.
    //! filters terms according to the provided pattern
    Pk_rsd(const Pk_one_loop& Pk, const GiNaC::symbol& mu_, const filter_list pt_, const GiNaC_symbol_set sy_);

    //! destructor is default
    ~Pk_rsd() = default;


    // INTERNAL API

  protected:

    //! filter a Pk_db into a destination Pk_rsd_group
    void filter(Pk_rsd_group& dest, const Pk_one_loop_impl::Pk_db& source);


    // ACCESSORS

  public:

    //! get tree-level terms
    const Pk_rsd_group& get_tree() const { return this->Ptree; }

    //! get 13 terms
    const Pk_rsd_group& get_13() const { return this->P13; }

    //! get 22 terms
    const Pk_rsd_group& get_22() const { return this->P22; }


    // SERVICES

  public:

    //! write self to stream
    void write(std::ostream& out) const;


    // INTERNAL DATA

  private:

    //! cache mu parameter
    const GiNaC::symbol mu;


    // FILTER/PATTERN

    //! cache pattern used to filter input expression
    const filter_list pattern;


    // DATABASES

    //! tree contributions
    Pk_rsd_group Ptree;

    //! 13 contributions
    Pk_rsd_group P13;

    //! 22 contributions
    Pk_rsd_group P22;

  };


//! stream insertion operator
std::ostream& operator<<(std::ostream& str, const Pk_rsd& obj);


#endif //LSSEFT_ANALYTIC_PK_RSD_H
