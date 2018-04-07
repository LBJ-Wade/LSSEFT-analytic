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

#include "Bk_cross_products.h"


namespace cross_product_impl
  {

    kernel_product relabel_kernel_product(const fourier_kernel_impl::kernel& ker1, GiNaC::symbol k1,
                                          const fourier_kernel_impl::kernel& ker2, GiNaC::symbol k2,
                                          const fourier_kernel_impl::kernel& ker3, GiNaC::symbol k3,
                                          const Wick::Wick_data& data, service_locator& loc)
      {
        // extract data from each kernel
        const auto& K1 = ker1.get_kernel();
        const auto& K2 = ker2.get_kernel();
        const auto& K3 = ker3.get_kernel();

        const auto& rm1 = ker1.get_substitution_list();
        const auto& rm2 = ker2.get_substitution_list();
        const auto& rm3 = ker3.get_substitution_list();

        const auto& loops = data.get_loop_momenta();

        // STEP 1. PERFORM RELABELLING

        // before taking the product K1*K2*K3 we must relabel indices and Rayleigh momenta

        // get substitutions produced by the Wick contractions
        const auto& subs_maps = data.get_substitution_rules();

        // there should be three substitution maps, one each for K1, K2, K3. If not, complain
        if(subs_maps.size() != 3)
          throw exception(ERROR_INCORRECT_SUBMAP_SIZE, exception_code::Bk_error);

        // STEP 1a. Merge (and relabel) Rayleigh momenta

        // first, generate an empty list of Rayleigh momenta. This will eventually hold the merged list.
        GiNaC::exmap Rayleigh_list;
        // second, build a list of labels that we must avoid when generating new momenta.
        // this should contain the external momentum 'k' and any loop momenta
        GiNaC_symbol_set reserved{k1, k2, k3};
        std::copy(loops.begin(), loops.end(), std::inserter(reserved, reserved.begin()));

        // merge the Rayleigh momenta from each kernel into Rayleigh_list, applying any necessary substitutions
        // and avoiding any labels in 'reserved'.
        // The result of the merge is a relabelling rule for each kernel.
        using Rayleigh::merge_Rayleigh_lists;
        auto Ray_remap1 = merge_Rayleigh_lists(rm1, Rayleigh_list, reserved, subs_maps[0], loc);
        auto Ray_remap2 = merge_Rayleigh_lists(rm2, Rayleigh_list, reserved, subs_maps[1], loc);
        auto Ray_remap3 = merge_Rayleigh_lists(rm3, Rayleigh_list, reserved, subs_maps[2], loc);

        // STEP 1b. Perform relabellings from Wick contraction and Rayleigh relabelling

        auto K1_remap = K1.subs(subs_maps[0]).subs(Ray_remap1);
        auto K2_remap = K2.subs(subs_maps[1]).subs(Ray_remap2);
        auto K3_remap = K3.subs(subs_maps[2]).subs(Ray_remap3);

        // STEP 2. Build the product K1*K2, relabelling any common indices in the process

        auto temp = relabel_index_product(K1_remap, K2_remap, loc);
        auto K = relabel_index_product(temp, K3_remap, loc);

        return std::make_pair(K, Rayleigh_list);
      }


    void simplify_kernel_product(kernel_product& product, GiNaC::symbol k1, GiNaC::symbol k2, GiNaC::symbol k3,
                                 service_locator& loc)
      {
        auto& K = product.first;
        auto& Rayleigh_list = product.second;

        // STEP 1. REMOVE TRIVIAL RAYLEIGH RULES

        // a trivial rule is a replacement such as p -> q. These cause no problems.
        // the mappings we need to worry about are those of the form p -> k-q

        using Rayleigh::remove_Rayleigh_trivial;
        auto Rayleigh_triv = remove_Rayleigh_trivial(Rayleigh_list);
        K = K.subs(Rayleigh_triv);


        // STEP 2. SIMPLIFY DOT PRODUCTS WHERE POSSIBLE

        GiNaC::scalar_products dotp;
        dotp.add(k1, k1, k1*k1);
        dotp.add(k2, k2, k2*k2);
        dotp.add(k3, k3, k3*k3);
        dotp.add(k1, k2, (k3*k3 - k1*k1 - k2*k2)/2);
        dotp.add(k1, k3, (k2*k2 - k1*k1 - k3*k3)/2);
        dotp.add(k2, k3, (k1*k1 - k2*k2 - k3*k3)/2);
        // k.l, l.l and other inner products are supposed to be picked up later by
        // loop integral transformations

        K = simplify_index(K, dotp, Rayleigh_list, loc);


        // STEP 3. PRUNE RAYLEIGH RULES TO REMOVE UNNEEDED ENTRIES

        // some mappings may now be redundant as a result of simplification of dot products
        using Rayleigh::prune_Rayleigh_list;
        prune_Rayleigh_list(Rayleigh_list, K);
      }

  }   // namespace cross_product_impl
