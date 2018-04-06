//
// Created by David Seery on 25/08/2017.
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

#ifndef LSSEFT_ANALYTIC_PK_ONE_LOOP_H
#define LSSEFT_ANALYTIC_PK_ONE_LOOP_H


#include <iostream>

#include "lib/fourier_kernel.h"
#include "lib/expression_databases/tree_db.h"
#include "lib/expression_databases/oneloop_db.h"
#include "lib/correlators/detail/Pk_cross_products.h"
#include "lib/correlators/detail/relabel_product.h"
#include "lib/correlators/detail/Rayleigh_momenta.h"

#include "services/symbol_factory.h"


//! Pk_one_loop understands how to construct the one-loop power spectrum from a set of Fourier kernels
class Pk_oneloop
  {

    // TYPES
    
  public:

    //! pull in tree_db as our database type for tree terms
    using Pk_tree_db = tree_db;

    //! pull in oneloop_db as our database type for 1-loop terms
    using Pk_oneloop_db = oneloop_db;


    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor accepts two Fourier kernels and the corresponding momentum labels
    template <unsigned int N1, unsigned int N2>
    Pk_oneloop(std::string n_, std::string t_, const fourier_kernel<N1>& ker1, const fourier_kernel<N2>& ker2,
                GiNaC::symbol k_, service_locator& lc_);

    //! copy constructor
    Pk_oneloop(const Pk_oneloop& obj);
    
    //! destructor is default
    ~Pk_oneloop() = default;
    
    
    // BUILD POWER SPECTRA EXPRESSIONS
    
  protected:

    //! build tree power spectrum
    template <typename Kernel1, typename Kernel2>
    void build_tree(const Kernel1& ker1, const Kernel2& ker2);
    
    //! build 13 power spectrum
    template <typename Kernel1, typename Kernel2>
    void build_13(const Kernel1& ker1, const Kernel2& ker2);
    
    //! build 22 power spectrum
    template <typename Kernel1, typename Kernel2>
    void build_22(const Kernel1& ker1, const Kernel2& ker2);

    //! perform reduction of angular integrals
    void perform_angular_reduction();

    //! clear all reduced angular integrals
    void clear_angular_reductions();


    // ADD OR SUBTRACT CORRELATION FUNCTIONS

  public:

    //! increment using 2nd correlation function
    Pk_oneloop& operator+=(const Pk_oneloop& obj);


    // EXTRACT EXPRESSIONS
    
  public:
    
    //! get tree power spectrum
    const Pk_tree_db& get_tree() const { return this->Ptree; }
    
    //! get 13 power spectrum
    const Pk_oneloop_db& get_13() const { return this->P13; }
    
    //! get 22 power spectrum
    const Pk_oneloop_db& get_22() const { return this->P22; }


    // TRANSFORMATIONS

  public:

    //! apply simplifications
    void simplify(const GiNaC::exmap& map);

    //! canonicalize external momenta
    void canonicalize_external_momenta();


    // SERVICES

  public:

    //! write self to stream
    void write(std::ostream& out) const;

    //! write Mathematica script for loop integrals
    void write_Mathematica(std::ostream& out) const;

    //! get tag
    const std::string& get_tag() const { return this->tag; }


    // FRIEND DECLARATIONS

    friend Pk_oneloop operator+(const Pk_oneloop& a, const Pk_oneloop& b);


    // INTERNAL DATA
    
  private:
    
    // SERVICES
    
    //! cache reference to service locator
    service_locator& loc;
    
    
    // RESERVED SYMBOLS

    //! cache momentum label k
    const GiNaC::symbol k;


    // METADATA

    //! cache name
    std::string name;

    //! cache tag
    std::string tag;


    // POWER SPECTRUM EXPRESSIONS
    
    //! expression for tree power spectrum
    Pk_tree_db Ptree;
    
    //! expression for 13 power spectrum
    Pk_oneloop_db P13;
    
    //! expression for 22 power spectrum
    Pk_oneloop_db P22;
  
  };


//! perform stream insertion
std::ostream& operator<<(std::ostream& str, const Pk_oneloop& obj);


//! add two power spectra
Pk_oneloop operator+(const Pk_oneloop& a, const Pk_oneloop& b);


template <unsigned int N1, unsigned int N2>
Pk_oneloop::Pk_oneloop(std::string n_, std::string t_, const fourier_kernel<N1>& ker1, const fourier_kernel<N2>& ker2,
                         GiNaC::symbol k_, service_locator& lc_)
  : name(std::move(n_)),
    tag(std::move(t_)),
    k(std::move(k_)),
    loc(lc_)
  {
    static_assert(N1 >= 3, "To construct a 1-loop power spectrum requires ker1 to be a Fourier kernel of third-order or above");
    static_assert(N2 >= 3, "To construct a 1-loop power spectrum requires ker2 to be a Fourier kernel of third-order or above");

    // build power spectrum components from products of ker1 and ker2
    this->build_tree(ker1, ker2);
    this->build_13(ker1, ker2);
    this->build_22(ker1, ker2);

    this->perform_angular_reduction();
  }


template <typename Kernel1, typename Kernel2>
void Pk_oneloop::build_tree(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db = ker1.order(1);
    const auto ker2_db = ker2.order(1);

    cross_product(ker1_db, ker2_db, this->k, this->Ptree, this->loc);
  }


template <typename Kernel1, typename Kernel2>
void Pk_oneloop::build_13(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db1 = ker1.order(1);
    const auto ker1_db3 = ker1.order(3);
    
    const auto ker2_db1 = ker2.order(1);
    const auto ker2_db3 = ker2.order(3);

    cross_product(ker1_db1, ker2_db3, this->k, this->P13, this->loc);
    cross_product(ker1_db3, ker2_db1, this->k, this->P13, this->loc);
  }


template <typename Kernel1, typename Kernel2>
void Pk_oneloop::build_22(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db2 = ker1.order(2);
    const auto ker2_db2 = ker2.order(2);

    cross_product(ker1_db2, ker2_db2, this->k, this->P22, this->loc);
  }


#endif //LSSEFT_ANALYTIC_PK_ONE_LOOP_H
