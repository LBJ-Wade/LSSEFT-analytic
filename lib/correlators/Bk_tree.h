//
// Created by David Seery on 05/04/2018.
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

#ifndef LSSEFT_ANALYTIC_BK_TREE_H
#define LSSEFT_ANALYTIC_BK_TREE_H


#include <iostream>

#include "lib/fourier_kernel.h"
#include "lib/expression_databases/tree_db.h"
#include "lib/correlators/detail/contractions.h"
#include "lib/correlators/detail/relabel_product.h"
#include "lib/correlators/detail/Rayleigh_momenta.h"

#include "services/symbol_factory.h"


//! Bk_tree understands how to construct the tree-level bispectrum from a set of Fourier kernels
class Bk_tree
  {

    // TYPES

  public:

    // pull in tree_db as our database type
    using Bk_db = tree_db;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor accepts three Fourier kernels and their corresponding momentum labels
    template <unsigned int N1, unsigned int N2, unsigned int N3>
    Bk_tree(std::string n_, std::string t_, const fourier_kernel<N1>& ker1, GiNaC::symbol k1_,
            const fourier_kernel<N2>& ker2, GiNaC::symbol k2_,
            const fourier_kernel<N3>& ker3, GiNaC::symbol k3_, service_locator& lc_);

    //! copy constructor
    Bk_tree(const Bk_tree& obj);

    //! destructor is default
    ~Bk_tree() = default;


    // BUILD BISPECTRUM EXPRESSIONS

    //! generic algorithm to construct cross-products of kernels
    //! (assumed to be of a single order, but the algorithm doesn't enforce that)
    template <typename Kernel1, typename Kernel2, typename Kernel3>
    void cross_product(const Kernel1& ker1, const Kernel2& ker2, const Kernel3& ker3, Bk_db& db);

    //! build tree-level bispectrum
    template <typename Kernel1, typename Kernel2, typename Kernel3>
    void build_tree(const Kernel1& ker1, const Kernel2& ker2, const Kernel3& ker3);


    // ADD OR SUBTRACT CORRELATION FUNCTIONS

  public:

    //! increment using 2nd correlation function
    Bk_tree& operator+=(const Bk_tree& obj);


    // EXTRACT EXPRESSIONS

  public:

    //! get tree power spectrum
    const Bk_db& get_tree() const { return this->Btree; }


    // SERVICES

  public:

    //! write self to stream
    void write(std::ostream& out) const;

    //! write Mathematica script
    void write_Mathematica(std::ostream& out) const;

    //! get tag
    const std::string& get_tag() const { return this->tag; }


    // FRIEND DECLARATIONS

    friend Bk_tree operator+(const Bk_tree& a, const Bk_tree& b);


    // INTERNAL DATA

  private:

    // SERVICES

    //! cache reference to service locator
    service_locator& loc;


    // RESERVED SYMBOLS

    //! cache momentum labels k1, k2, k3
    //! note that we keep track of the ordering, because we don't know that the fields appearing in the
    //! correlation function are permutable
    const GiNaC::symbol k1;
    const GiNaC::symbol k2;
    const GiNaC::symbol k3;


    // METADATA

    //! cache name
    std::string name;

    //! cache tag
    std::string tag;


    // BISPECTRUM EXPRESSIONS

    //! expressions for tree-level bispectrum
    Bk_db Btree;

  };


//! perform stream insertion
std::ostream& operator<<(std::ostream& str, const Bk_tree& obj);


//! add two bispectra
Bk_tree operator+(const Bk_tree& a, const Bk_tree& b);


template <unsigned int N1, unsigned int N2, unsigned int N3>
Bk_tree::Bk_tree(std::string n_, std::string t_, const fourier_kernel<N1>& ker1, GiNaC::symbol k1_,
                 const fourier_kernel<N2>& ker2, GiNaC::symbol k2_, const fourier_kernel<N3>& ker3, GiNaC::symbol k3_,
                 service_locator& lc_)
  : name(std::move(n_)),
    tag(std::move(t_)),
    k1(k1_),
    k2(k2_),
    k3(k3_),
    loc(lc_)
  {
    static_assert(N1 >= 1, "To construct a tree-level bispectrum requires ker1 to be Fourier kernel of first-order or above");
    static_assert(N2 >= 1, "To construct a tree-level bispectrum requires ker2 to be Fourier kernel of first-order or above");
    static_assert(N3 >= 1, "To construct a tree-level bispectrum requires ker3 to be Fourier kernel of first-order or above");

    // build bispectrum components from products of ker1, ker2 and ker3
    this->build_tree(ker1, ker2, ker3);
  }


template <typename Kernel1, typename Kernel2, typename Kernel3>
void Bk_tree::cross_product(const Kernel1& ker1, const Kernel2& ker2, const Kernel3& ker3, Bk_db& db)
  {
    // multiply out all terms in ker1, ker2 and ker3, using the insertion operator 'ins'
    // to store the results in a suitable Bk database
  }


template <typename Kernel1, typename Kernel2, typename Kernel3>
void Bk_tree::build_tree(const Kernel1& ker1, const Kernel2& ker2, const Kernel3& ker3)
  {
    const auto ker1_db1 = ker1.order(1);
    const auto ker1_db2 = ker1.order(2);

    const auto ker2_db1 = ker2.order(1);
    const auto ker2_db2 = ker2.order(2);

    const auto ker3_db1 = ker3.order(1);
    const auto ker3_db2 = ker3.order(2);

    this->cross_product(ker1_db1, ker2_db1, ker3_db2, this->Btree);
    this->cross_product(ker1_db1, ker2_db2, ker3_db1, this->Btree);
    this->cross_product(ker1_db2, ker2_db1, ker3_db1, this->Btree);
  }


#endif //LSSEFT_ANALYTIC_BK_TREE_H
