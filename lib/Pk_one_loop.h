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
#include <memory>

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

      public:

        //! a loop_pair is a doublet that groups a raw integral with its reduced form
        using loop_pair = std::pair< std::unique_ptr<loop_integral>, std::unique_ptr<one_loop_reduced_integral> >;

      protected:

        //! database is a set of loop_pairs, keyed by loop_integral_key
        using db_type = std::unordered_map< loop_integral_key, loop_pair >;

      public:

        //! iterator
        using iterator = db_type::iterator;
        using const_iterator = db_type::const_iterator;


        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor is default
        Pk_db() = default;

        //! copy constructor
        Pk_db(const Pk_db& obj);

        //! destructor is default
        ~Pk_db() = default;


        // ITERATORS

      public:

        //! iterators
        iterator       begin()        { return this->db.begin(); }
        iterator       end()          { return this->db.end(); }
        const_iterator begin()  const { return this->db.cbegin(); }
        const_iterator end()    const { return this->db.cend(); }

        const_iterator cbegin() const { return this->db.cbegin(); }
        const_iterator cend()   const { return this->db.cend(); }


        // EMPLACE AN ELEMENT

      public:

        //! emplace an element
        void emplace(std::unique_ptr<loop_integral> elt);


        // FLUSH REDUCED INTEGRALS

      public:

        //! clear all reduced integrals
        void clear_reduced_integrals();


        // ARITHMETIC

      public:

        //! increment
        Pk_db& operator+=(const Pk_db& obj);


        // TRANSFORMATIONS

      public:

        //! reduce angular integrals
        void reduce_angular_integrals(service_locator& loc, bool symmetrize);

        //! apply simplification map
        void simplify(const GiNaC::exmap& map);

        //! canonicalize external momenta
        void canonicalize_external_momenta();

        //! prune empty records
        void prune();


        // SERVICES

      public:

        //! write self to steam
        void write(std::ostream& out) const;

        //! write Mathematica-format expression
        void write_Mathematica(std::ostream& out, std::string symbol, bool do_dx) const;

        //! compute UV limit
        GiNaC::ex get_UV_limit(unsigned int order=2) const;


        // INTERNAL DATA

      protected:

        //! database of loop integrals
        db_type db;

      };

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
    Pk_one_loop(std::string n_, std::string t_, const fourier_kernel<N1>& ker1, const fourier_kernel<N2>& ker2,
                GiNaC::symbol k_, service_locator& lc_);

    //! copy constructor
    Pk_one_loop(const Pk_one_loop& obj);
    
    //! destructor is default
    ~Pk_one_loop() = default;
    
    
    // BUILD POWER SPECTRA EXPRESSIONS
    
  protected:
    
    //! generic algorithm to construct cross-products of kernels
    //! (assumed to be of a single order, but the algorithm doesn't enforce that)
    template <typename Kernel1, typename Kernel2>
    void cross_product(const Kernel1& ker1, const Kernel2& ker2, Pk_db& db);
    
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
    Pk_one_loop& operator+=(const Pk_one_loop& obj);


    // EXTRACT EXPRESSIONS
    
  public:
    
    //! get tree power spectrum
    const Pk_db& get_tree() const { return this->Ptree; }
    
    //! get 13 power spectrum
    const Pk_db& get_13() const { return this->P13; }
    
    //! get 22 power spectrum
    const Pk_db& get_22() const { return this->P22; }


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


    // INTERNAL DATA
    
  private:
    
    // SERVICES
    
    //! cache reference to symbol factory
    service_locator& loc;
    
    
    // RESERVED SYMBOLS

    //! cache momentum label k
    const GiNaC::symbol k;


    // METADATA

    //! cache name
    const std::string name;

    //! cache tag
    const std::string tag;


    // POWER SPECTRUM EXPRESSIONS
    
    //! expression for tree power spectrum
    Pk_db Ptree;
    
    //! expression for 13 power spectrum
    Pk_db P13;
    
    //! expression for 22 power spectrum
    Pk_db P22;
  
  };


//! perform stream insertion
std::ostream& operator<<(std::ostream& str, const Pk_one_loop& obj);

//! add two power spectrum
inline Pk_one_loop operator+(const Pk_one_loop& a, const Pk_one_loop& b)
  {
    Pk_one_loop c{a};
    c += b;

    return a;
  }


template <unsigned int N1, unsigned int N2>
Pk_one_loop::Pk_one_loop(std::string n_, std::string t_, const fourier_kernel<N1>& ker1, const fourier_kernel<N2>& ker2,
                         GiNaC::symbol k_, service_locator& lc_)
  : name(std::move(n_)),
    tag(std::move(t_)),
    k(std::move(k_)),
    loc(lc_)
  {
    static_assert(N1 >= 3, "To construct a one-loop power spectrum requires a Fourier kernel of third-order or above");
    static_assert(N2 >= 3, "To construct a one-loop power spectrum requires a Fourier kernel of third-order or above");

    // build power spectrum components from products of ker1 and ker2
    this->build_tree(ker1, ker2);
    this->build_13(ker1, ker2);
    this->build_22(ker1, ker2);

    this->perform_angular_reduction();
  }


template <typename Kernel1, typename Kernel2>
void Pk_one_loop::cross_product(const Kernel1& ker1, const Kernel2& ker2, Pk_db& db)
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
                                      detail::contractions::kext_group<2>{ this->k, -this->k }, this->loc);
        
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
                auto Ray_remap1 = merge_Rayleigh_lists(rm1, Rayleigh_list, reserved, subs_maps[0], this->loc);
                auto Ray_remap2 = merge_Rayleigh_lists(rm2, Rayleigh_list, reserved, subs_maps[1], this->loc);
                
                // perform all relabellings
                auto K1_remap = K1.subs(subs_maps[0]).subs(Ray_remap1);
                auto K2_remap = K2.subs(subs_maps[1]).subs(Ray_remap2);
                
                // relabel indices
                using detail::relabel_index_product;
                auto K = relabel_index_product(K1_remap, K2_remap, this->loc);
                
                using detail::remove_Rayleigh_trivial;
                auto Rayleigh_triv = remove_Rayleigh_trivial(Rayleigh_list);
                K = K.subs(Rayleigh_triv);
                
                // simplify dot products where possible
                GiNaC::scalar_products dotp;
                dotp.add(this->k, this->k, this->k*this->k);
                // k.l, l.l and other inner products are supposed to be picked up later by
                // loop integral transformations

                K = simplify_index(K, dotp, Rayleigh_list, this->loc);

                // prune Rayleigh list to remove momenta that have dropped out
                using detail::prune_Rayleigh_list;
                prune_Rayleigh_list(Rayleigh_list, K);
                
                if(static_cast<bool>(K != 0))
                  {
                    auto elt =
                      std::make_unique<loop_integral>(tm1*tm2, K, data.get_Wick_string(), loops,
                                                      GiNaC_symbol_set{this->k}, Rayleigh_list, this->loc);
                    db.emplace(std::move(elt));
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

    this->cross_product(ker1_db, ker2_db, this->Ptree);
  }


template <typename Kernel1, typename Kernel2>
void Pk_one_loop::build_13(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db1 = ker1.order(1);
    const auto ker1_db3 = ker1.order(3);
    
    const auto ker2_db1 = ker2.order(1);
    const auto ker2_db3 = ker2.order(3);

    this->cross_product(ker1_db1, ker2_db3, this->P13);
    this->cross_product(ker1_db3, ker2_db1, this->P13);
  }


template <typename Kernel1, typename Kernel2>
void Pk_one_loop::build_22(const Kernel1& ker1, const Kernel2& ker2)
  {
    const auto ker1_db2 = ker1.order(2);
    const auto ker2_db2 = ker2.order(2);

    this->cross_product(ker1_db2, ker2_db2, this->P22);
  }


#endif //LSSEFT_ANALYTIC_PK_ONE_LOOP_H
