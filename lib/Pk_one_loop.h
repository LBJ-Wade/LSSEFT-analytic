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
#include <unordered_map>

#include "fourier_kernel.h"
#include "loop_integral.h"
#include "one_loop_reduced_integral.h"
#include "detail/contractions.h"
#include "detail/relabel_product.h"
#include "detail/Rayleigh_momenta.h"

#include "services/symbol_factory.h"


namespace Pk_one_loop_impl
  {

    //! a Pk_db is a container to the loop integrals generated as part of a power spectrum computation
    class Pk_db
      {

        // TYPES

      protected:

        //! database is a set of pairs that link the raw loop integrals with their
        //! reduced forms
        using db_type = std::vector< std::pair< std::unique_ptr<loop_integral>,
                                                std::unique_ptr<one_loop_reduced_integral> > >;


        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor is default
        Pk_db() = default;

        //! destructor is default
        ~Pk_db() = default;


        // EMPLACE AN ELEMENT

      public:

        //! emplace an element
        template <typename... Args>
        void emplace(Args&&... args);


        // TRANSFORMATIONS

      public:

        //! reduce angular integrals
        void reduce_angular_integrals();


        // SERVICES

      public:

        //! write self to steam
        void write(std::ostream& out) const;


        // INTERNAL DATA

      protected:

        //! database of loop integrals
        db_type db;

      };


    template <typename... Args>
    void Pk_db::emplace(Args&& ... args)
      {
        this->db.emplace_back( std::make_pair(std::forward<Args>(args)..., std::unique_ptr<one_loop_reduced_integral>()) );
      }

  }   // namespace Pk_one_loop_impl


//! perform stream insertion
std::ostream& operator<<(std::ostream& str, const Pk_one_loop_impl::Pk_db& obj);


//! Pk_one_loop understands how to construct the one-loop power spectrum from a set of Fourier kernels
class Pk_one_loop
  {

    // TYPES
    
  public:

    // pull in Pk_db
    using Pk_db = Pk_one_loop_impl::Pk_db;


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
    
    //! generic algorithm to construct cross-products of kernels
    //! (assumed to be of a single order, but the algorithm doesn't enforce that)
    template <typename Kernel1, typename Kernel2, typename InsertOperator>
    void cross_product(const Kernel1& ker1, const Kernel2& ker2, InsertOperator ins);
    
    //! build tree power spectrum
    template <typename Kernel1, typename Kernel2>
    void build_tree(const Kernel1& ker1, const Kernel2& ker2);
    
    //! build 13 power spectrum
    template <typename Kernel1, typename Kernel2>
    void build_13(const Kernel1& ker1, const Kernel2& ker2);
    
    //! build 22 power spectrum
    template <typename Kernel1, typename Kernel2>
    void build_22(const Kernel1& ker1, const Kernel2& ker2);


    // EXTRACT EXPRESSIONS
    
  public:
    
    //! get tree power spectrum
    const Pk_db& get_tree() const { return this->Ptree; }
    
    //! get 13 power spectrum
    const Pk_db& get_13() const { return this->P13; }
    
    //! get 22 power spectrum
    const Pk_db& get_22() const { return this->P22; }
    
    
    // INTERNAL DATA
    
  private:
    
    // SERVICES
    
    //! cache reference to symbol factory
    symbol_factory& sf;
    
    
    // RESERVED SYMBOLS

    //! cache momentum label k
    GiNaC::symbol k;


    // POWER SPECTRUM EXPRESSIONS
    
    //! expression for tree power spectrum
    Pk_db Ptree;
    
    //! expression for 13 power spectrum
    Pk_db P13;
    
    //! expression for 22 power spectrum
    Pk_db P22;
  
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

    this->P13.reduce_angular_integrals();
    this->P22.reduce_angular_integrals();
  }


template <typename Kernel1, typename Kernel2, typename InsertOperator>
void Pk_one_loop::cross_product(const Kernel1& ker1, const Kernel2& ker2, InsertOperator ins)
  {
    // multiply out all terms in ker1 and ker2, using the insertion operator 'ins'
    // to store the results in a suitable Pk database
    
    for(auto t1 = ker1.cbegin(); t1 != ker1.cend(); ++t1)
      {
        for(auto t2 = ker2.cbegin(); t2 != ker2.cend(); ++t2)
          {
            const auto& tm1 = t1->second->get_time_function();
            const auto& tm2 = t2->second->get_time_function();
        
            const auto& K1 = t1->second->get_kernel();
            const auto& K2 = t2->second->get_kernel();
        
            const auto& iv1 = t1->second->get_initial_value_set();
            const auto& iv2 = t2->second->get_initial_value_set();
            
            const auto& rm1 = t1->second->get_substitution_list();
            const auto& rm2 = t2->second->get_substitution_list();
        
            detail::contractions ctrs(detail::contractions::iv_group<2>{ iv1, iv2 },
                                      detail::contractions::kext_group<2>{ this->k, -this->k }, this->sf);
        
            const auto& Wicks = ctrs.get();
            for(const auto& W : Wicks)
              {
                const auto& data = *W;
                const auto& loops = data.get_loop_momenta();
                
                if(loops.size() > 1)
                  throw exception(ERROR_EXPECTED_ONE_LOOP_RESULT, exception_code::Pk_error);
    
                // before taking the product K1*K2 we must relabel indices in K2 if they clash with
                // K1, otherwise we will get nonsensical results
                const auto& subs_maps = data.get_substitution_rules();
                if(subs_maps.size() != 2)
                  throw exception(ERROR_INCORRECT_SUBMAP_SIZE, exception_code::Pk_error);
    
                // merge lists of Rayleigh rules together
                GiNaC::exmap Rayleigh_list;
                GiNaC_symbol_set reserved{k};
                std::copy(loops.begin(), loops.end(), std::inserter(reserved, reserved.begin()));
                
                using detail::merge_Rayleigh_lists;
                auto Ray_remap1 = merge_Rayleigh_lists(rm1, Rayleigh_list, reserved, subs_maps[0], this->sf);
                auto Ray_remap2 = merge_Rayleigh_lists(rm2, Rayleigh_list, reserved, subs_maps[1], this->sf);
                
                // perform all relabellings
                auto K1_remap = K1.subs(subs_maps[0]).subs(Ray_remap1);
                auto K2_remap = K2.subs(subs_maps[1]).subs(Ray_remap2);
                
                // relabel indices
                using detail::relabel_index_product;
                auto K = relabel_index_product(K1_remap, K2_remap, this->sf);
                
                using detail::remove_Rayleigh_trivial;
                auto Rayleigh_triv = remove_Rayleigh_trivial(Rayleigh_list);
                K = K.subs(Rayleigh_triv);
                
                // simplify dot products where possible
                GiNaC::scalar_products dotp;
                dotp.add(this->k, this->k, this->k*this->k);
                // k.l, l.l and other inner products are supposed to be picked up later by
                // loop integral transformations

                K = simplify_index(K, dotp);

                // prune Rayleigh list to remove momenta that have dropped out
                using detail::prune_Rayleigh_list;
                prune_Rayleigh_list(Rayleigh_list, K);
                
                if(static_cast<bool>(K != 0))
                  {
                    ins(tm1 * tm2, K, data.get_Wick_string(), loops, Rayleigh_list);
                  }
              }
          }
      }
  }


template <typename Kernel1, typename Kernel2>
void Pk_one_loop::build_tree(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db = ker1.order(1);
    const auto ker2_db = ker2.order(1);

    auto ins = [&](time_function t, GiNaC::ex K, GiNaC::ex ws, GiNaC_symbol_set lm, subs_list rm) -> void
      {
        this->Ptree.emplace(
          std::make_unique<loop_integral>(
            std::move(t), std::move(K), std::move(ws), std::move(lm), GiNaC_symbol_set{this->k}, std::move(rm), this->sf
          )
        );
      };
    
    this->cross_product(ker1_db, ker2_db, ins);
  }


template <typename Kernel1, typename Kernel2>
void Pk_one_loop::build_13(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db1 = ker1.order(1);
    const auto ker1_db3 = ker1.order(3);
    
    const auto ker2_db1 = ker2.order(1);
    const auto ker2_db3 = ker2.order(3);

    auto ins = [&](time_function t, GiNaC::ex K, GiNaC::ex ws, GiNaC_symbol_set lm, subs_list rm) -> void
      {
        this->P13.emplace(
          std::make_unique<loop_integral>(
            std::move(t), std::move(K), std::move(ws), std::move(lm), GiNaC_symbol_set{this->k}, std::move(rm), this->sf
          )
        );
      };
    
    this->cross_product(ker1_db1, ker2_db3, ins);
    this->cross_product(ker1_db3, ker2_db1, ins);
  }


template <typename Kernel1, typename Kernel2>
void Pk_one_loop::build_22(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db2 = ker1.order(2);
    const auto ker2_db2 = ker2.order(2);

    auto ins = [&](time_function t, GiNaC::ex K, GiNaC::ex ws, GiNaC_symbol_set lm, subs_list rm) -> void
      {
        this->P22.emplace(
          std::make_unique<loop_integral>(
            std::move(t), std::move(K), std::move(ws), std::move(lm), GiNaC_symbol_set{this->k}, std::move(rm), this->sf
          )
        );
      };
    
    this->cross_product(ker1_db2, ker2_db2, ins);
  }


#endif //LSSEFT_ANALYTIC_PK_ONE_LOOP_H
