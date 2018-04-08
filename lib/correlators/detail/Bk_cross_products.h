//
// Created by David Seery on 06/04/2018.
// --@@
// Copyright (c) 2018 University of Sussex. All rights reserved.
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
#ifndef LSSEFT_ANALYTIC_BK_CROSS_PRODUCTS_H
#define LSSEFT_ANALYTIC_BK_CROSS_PRODUCTS_H


#include <string>
#include <sstream>

#include "cross_products.h"

#include "lib/fourier_kernel.h"
#include "lib/correlators/detail/contractions.h"
#include "lib/correlators/detail/Rayleigh_momenta.h"
#include "lib/expression_databases/tree_db.h"
#include "lib/expression_databases/oneloop_db.h"

#include "services/service_locator.h"

#include "localizations/messages.h"
#include "shared/exceptions.h"


namespace cross_product_impl
  {

    // generic algorithm to perform cross products
    template <unsigned int N1, unsigned int N2, unsigned int N3, typename ProductHandler>
    void cross_product(const fourier_kernel<N1>& ker1, const fourier_kernel<N2>& ker2, const fourier_kernel<N3>& ker3,
                       unsigned int product_order, std::string error_msg, ProductHandler prod)
      {
        for(auto t1 = ker1.cbegin(); t1 != ker1.cend(); ++t1)
          {
            for(auto t2 = ker2.cbegin(); t2 != ker2.cend(); ++t2)
              {
                for(auto t3 = ker3.cbegin(); t3 != ker3.cend(); ++t3)
                  {
                    // check that these terms can be multiplied to product a result at the required level
                    if(t1->second->order() + t2->second->order() + t3->second->order() != product_order)
                      {
                        std::ostringstream msg;
                        msg << error_msg << " " << t1->second->order() << ", " << t2->second->order() << ", " << t3->second->order();
                        throw exception(msg.str(), exception_code::Pk_error);
                      }

                    // apply supplied product function
                    prod(*t1->second, *t2->second, *t3->second);
                  }
              }
          }
      }


    // construct product of three kernels, performing any necessary relabellings
    // returns a pair consisting of the relabelled product, and a set of Rayleigh rules
    kernel_product
    relabel_kernel_product(const fourier_kernel_impl::kernel& ker1, GiNaC::symbol k1,
                           const fourier_kernel_impl::kernel& ker2, GiNaC::symbol k2,
                           const fourier_kernel_impl::kernel& ker3, GiNaC::symbol k3,
                           const Wick::Wick_data& data, service_locator& loc);


    // perform in-place simplification of a kernel product
    void simplify_kernel_product(kernel_product& product, GiNaC::symbol k1, GiNaC::symbol k2, GiNaC::symbol k3,
                                 service_locator& loc);



    class Bk_tree_handler
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor captures database
        Bk_tree_handler(tree_db& db_, GiNaC::symbol k1_, GiNaC::symbol k2_, GiNaC::symbol k3_, service_locator& loc_)
          : db(db_),
            k1(std::move(k1_)),
            k2(std::move(k2_)),
            k3(std::move(k3_)),
            loc(loc_)
          {
          }

        //! destructor is default
        ~Bk_tree_handler() = default;

        // APPLY

      public:

        //! apply handler
        void operator()(const fourier_kernel_impl::kernel& ker1, const fourier_kernel_impl::kernel& ker2,
                        const fourier_kernel_impl::kernel& ker3)
          {
            // extract data from each kernel
            const auto& tm1 = ker1.get_time_function();
            const auto& tm2 = ker2.get_time_function();
            const auto& tm3 = ker3.get_time_function();

            const auto& iv1 = ker1.get_initial_value_set();
            const auto& iv2 = ker2.get_initial_value_set();
            const auto& iv3 = ker3.get_initial_value_set();

            // product the full set of (connected) Wick contractions from this pair of kernels
            Wick::contractions ctrs(Wick::contractions::iv_group<3>{ iv1, iv2, iv3 },
                                    Wick::contractions::kext_group<3>{ this->k1, this->k2, this->k3 }, this->loc);

            // loop through all Wick contractions, arranging the correct replacements and pushing
            // an appropriate product into the output database
            const auto& Wicks = ctrs.get();

            for(const auto& W : Wicks)
              {
                const auto& data = *W;
                const auto& loops = data.get_loop_momenta();

                if(loops.size() != 0)
                  {
                    std::ostringstream msg;
                    msg << ERROR_EXPECTED_TREE_RESULT << " " << loops.size();
                    throw exception(msg.str(), exception_code::Pk_error);
                  }

                // STEP 1. CONSTRUCT THE PRODUCT OF KERNELS
                auto product = relabel_kernel_product(ker1, this->k1, ker2, this->k2, ker3, this->k3, data, this->loc);

                // STEP 2. SIMPLIFY RESULT
                // (performs in-place modification of 'product' if needed)
                simplify_kernel_product(product, this->k1, this->k2, this->k3, this->loc);

                // STEP 3. PERFORM INSERTION
                auto& K = product.first;
                auto& Rayleigh_list = product.second;

                // check that no Rayleigh rules have been produced
                if(product.second.size() != 0)
                  {
                    const auto& map = product.second;
                    for(const auto& item : map)
                      {
                        std::cerr << item.first << " -> " << item.second << '\n';
                      }
                    std::cerr << "K = " << product.first << '\n';
                    throw exception(ERROR_EXPECTED_EMPTY_RAYLEIGH_LIST, exception_code::Pk_error);
                  }

                // if resulting kernel is not identically zero, push it to the database
                if(static_cast<bool>(K != 0))
                  {
                    auto elt =
                      std::make_unique<tree_expression>(tm1*tm2*tm3, K, data.get_Wick_string(),
                                                        GiNaC_symbol_set{this->k1, this->k2, this->k3}, this->loc);
                    this->db.emplace(std::move(elt));
                  }
              }
          }


        // INTERNAL DATA

      private:

        //! capture reference to service locator
        service_locator& loc;

        //! capture reference to database
        tree_db& db;

        //! momentum label k1
        GiNaC::symbol k1;

        //! momentum label k2
        GiNaC::symbol k2;

        //! momentum label k3
        GiNaC::symbol k3;

      };

  }   // namespace cross_product_impl


template <unsigned int N1, unsigned int N2, unsigned int N3>
void cross_product(const fourier_kernel<N1>& ker1, GiNaC::symbol k1,
                   const fourier_kernel<N2>& ker2, GiNaC::symbol k2,
                   const fourier_kernel<N3>& ker3, GiNaC::symbol k3,
                   tree_db& db, service_locator& loc)
  {
    using cross_product_impl::cross_product;
    using cross_product_impl::Bk_tree_handler;

    cross_product(ker1, ker2, ker3, 4, ERROR_EXPECTED_BK_TREE_PRODUCT, Bk_tree_handler(db, k1, k2, k3, loc));
  };

#endif //LSSEFT_ANALYTIC_BK_CROSS_PRODUCTS_H
