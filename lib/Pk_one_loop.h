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


#include <unordered_map>

#include "fourier_kernel.h"
#include "detail/contractions.h"

#include "services/symbol_factory.h"


//! Pk_one_loop understands how to construct the one-loop power spectrum from a set of
//! Fourier kernels
class Pk_one_loop
  {

    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor accepts two Fourier kernels and the corresponding momentum labels
    template <unsigned int N1, unsigned int N2>
    Pk_one_loop(const fourier_kernel<N1>& ker1, const fourier_kernel<N2>& ker2,
                const GiNaC::symbol& k_, symbol_factory& sf_);
    
    //! destructor is default
    ~Pk_one_loop() = default;
    
    
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

    
    // EXTRACT RAW EXPRESSIONS
    
  public:
    
    //! get tree power spectrum
    const GiNaC::ex& get_tree() const { return this->Ptree; }
    
    //! get 13 power spectrum
    const GiNaC::ex& get_13() const { return this->P13; }
    
    //! get 22 power spectrum
    const GiNaC::ex& get_22() const { return this->P22; }
    
    
    // INTERNAL DATA
    
  private:
    
    // SERVICES
    
    //! cache reference to symbol factory
    symbol_factory& sf;
    
    
    // POWER SPECTRUM COMPONENTS

    //! cache momentum label k
    GiNaC::symbol k;
    
    
    // POWER SPECTRUM EXPRESSIONS
    
    //! expression for tree power spectrum
    GiNaC::ex Ptree;
    
    //! expression for 13 power spectrum
    GiNaC::ex P13;
    
    //! expression for 22 power spectrum
    GiNaC::ex P22;
  
  };


template <unsigned int N1, unsigned int N2>
Pk_one_loop::Pk_one_loop(const fourier_kernel<N1>& ker1, const fourier_kernel<N2>& ker2,
                         const GiNaC::symbol& k_, symbol_factory& sf_)
  : k(k_),
    sf(sf_)
  {
    static_assert(N1 >= 3, "To construct a one-loop power spectrum requires a Fourier kernel of third-order or above");
    static_assert(N2 >= 3, "To construct a one-loop power spectrum requires a Fourier kernel of third-order or above");

    this->build_tree(ker1, ker2);
    this->build_13(ker1, ker2);
    this->build_22(ker1, ker2);
  }


template <typename Kernel1, typename Kernel2>
void Pk_one_loop::build_tree(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db = ker1.order(1);
    const auto ker2_db = ker2.order(1);
    
    // cross-multiply all terms in ker1_db and ker2_db
    for(auto t1 = ker1_db.cbegin(); t1 != ker1_db.cend(); ++t1)
      {
        for(auto t2 = ker2_db.cbegin(); t2 != ker2_db.cend(); ++t2)
          {
            const auto& tm1 = t1->second->get_time_function();
            const auto& tm2 = t2->second->get_time_function();
            
            const auto& K1 = t1->second->get_kernel();
            const auto& K2 = t2->second->get_kernel();
            
            const auto& iv1 = t1->second->get_initial_value_set();
            const auto& iv2 = t2->second->get_initial_value_set();

            detail::contractions ctrs(std::array<initial_value_set, 2>{iv1, iv2},
                                      std::array<GiNaC::ex, 2>{this->k, -this->k});

            const auto& Wicks = ctrs.get();

            for(const auto& W : Wicks)
              {
                const auto& data = *W;

                this->Ptree += tm1*tm2 * K1*K2 * data.get_Wick_string();
              }
          }
      }
  }


template <typename Kernel1, typename Kernel2>
void Pk_one_loop::build_13(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db1 = ker1.order(1);
    const auto ker1_db3 = ker1.order(3);
    
    const auto ker2_db1 = ker2.order(1);
    const auto ker2_db3 = ker2.order(3);
  }


template <typename Kernel1, typename Kernel2>
void Pk_one_loop::build_22(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db2 = ker1.order(2);
    const auto ker2_db2 = ker2.order(2);
  };


#endif //LSSEFT_ANALYTIC_PK_ONE_LOOP_H
