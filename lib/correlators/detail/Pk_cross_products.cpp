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

#include "Pk_cross_products.h"


namespace cross_product_impl
  {

    // construct product of two kernels, performing any necessary relabellings
    // returns a pair consisting of the relabelled product, and a set of Rayleigh rules
    std::pair< GiNaC::ex, GiNaC::exmap >
    relabel_kernel_product(const fourier_kernel_impl::kernel& ker1, const fourier_kernel_impl::kernel& ker2,
                           GiNaC::symbol k, const Wick::Wick_data& data, service_locator& loc)
      {
        // extract data from each kernel
        const auto& tm1 = ker1.get_time_function();
        const auto& tm2 = ker2.get_time_function();

        const auto& K1 = ker1.get_kernel();
        const auto& K2 = ker2.get_kernel();

        const auto& rm1 = ker1.get_substitution_list();
        const auto& rm2 = ker2.get_substitution_list();

        const auto& loops = data.get_loop_momenta();

        // STEP 1. PERFORM RELABELLING

        // before taking the product K1*K2 we must relabel indices and Rayleigh momenta in K2 if they clash
        // with K1, otherwise we will get nonsensical results

        // get substitutions produced by the Wick contractions
        const auto& subs_maps = data.get_substitution_rules();

        // there should be two substitution maps, one for K1 and one for K2. If not, complain
        if(subs_maps.size() != 2)
          throw exception(ERROR_INCORRECT_SUBMAP_SIZE, exception_code::Pk_error);

        // STEP 1a. Merge (and relabel) Rayleigh momenta

        // first, generate an empty list of Rayleigh momenta. This will eventually hold the merged list.
        GiNaC::exmap Rayleigh_list;
        // second, build a list of labels that we must avoid when generating new momenta.
        // this should contain the external momentum 'k' and any loop momenta
        GiNaC_symbol_set reserved{k};
        std::copy(loops.begin(), loops.end(), std::inserter(reserved, reserved.begin()));

        // merge the Rayleigh momenta from each kernel into Rayleigh_list, applying any necessary substitutions
        // and avoiding any labels in 'reserved'.
        // The result of the merge is a relabelling rule for each kernel.
        using Rayleigh::merge_Rayleigh_lists;
        auto Ray_remap1 = merge_Rayleigh_lists(rm1, Rayleigh_list, reserved, subs_maps[0], loc);
        auto Ray_remap2 = merge_Rayleigh_lists(rm2, Rayleigh_list, reserved, subs_maps[1], loc);

        // STEP 1b. Perform relabellings from Wick contraction and Rayleigh relabelling

        auto K1_remap = K1.subs(subs_maps[0]).subs(Ray_remap1);
        auto K2_remap = K2.subs(subs_maps[1]).subs(Ray_remap2);

        // STEP 2. Build the product K1*K2, relabelling any common indices in the process

        auto K = relabel_index_product(K1_remap, K2_remap, loc);

        return std::make_pair(K, Rayleigh_list);
      }


    void simplify_kernel_product(kernel_product& product, GiNaC::symbol k, service_locator& loc)
      {
        auto& K = product.first;
        auto& Rayleigh_list = product.second;

        // STEP 1. REMOVE TRIVIAL RAYLEIGH RULES

        // a trivial rule is a replacement such as p -> q. These cause no problems.
        // the mappings we need to worry about are those of the form p -> k-q

        using Rayleigh::remove_Rayleigh_trivial;
        auto Rayleigh_triv = remove_Rayleigh_trivial(Rayleigh_list);
        K = K.subs(Rayleigh_triv);


        // STEP 3. SIMPLIFY DOT PRODUCTS WHERE POSSIBLE

        // simplify k.k to k^2
        GiNaC::scalar_products dotp;
        dotp.add(k, k, k*k);
        // k.l, l.l and other inner products are supposed to be picked up later by
        // loop integral transformations

        K = simplify_index(K, dotp, Rayleigh_list, loc);


        // STEP 4. PRUNE RAYLEIGH RULES TO REMOVE UNNEEDED ENTRIES

        // some mappings may now be redundant as a result of simplification of dot products

        using Rayleigh::prune_Rayleigh_list;
        prune_Rayleigh_list(Rayleigh_list, K);
      }

  }   // namespace cross_product_impl
