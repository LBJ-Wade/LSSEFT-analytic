//
// Created by David Seery on 29/08/2017.
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

#include "relabel_product.h"


namespace detail
  {

    GiNaC::ex relabel_index_product(const GiNaC::ex& a, const GiNaC::ex& b, symbol_factory& sf)
      {
        GiNaC::exmap idx_map;
        
        // need only relabel indices that occur >= 2 times; can allow single occurrences to be contracted
        const auto& a_idxs = get_expr_indices(a, 2);
        const auto& b_idxs = get_expr_indices(b, 2);
        
        for(const auto& idx : b_idxs)
        {
          // does this index already exist in the a set? if not, nothing to do
          if(a_idxs.find(idx) == a_idxs.end()) continue;
          
          // otherwise, manufacture a new index
          // note that we remap the index *value* (ie. its symbolic label), not the index itself
            const auto relabel = sf.make_unique_index();
            idx_map[idx] = GiNaC::ex_to<GiNaC::symbol>(relabel.get_value());
        }
        
        return a * b.subs(idx_map);
      }
    
    
    GiNaC::exmap
    merge_Rayleigh_rules(GiNaC::exmap& dest, const GiNaC::exmap& source, const GiNaC_symbol_set& reserved,
                         const GiNaC::exmap& subs_rules, symbol_factory& sf)
      {
        GiNaC::exmap mma_map;
        
        // check whether any of the substitution rules in 'source' already exist in 'dest'
        // if they do, add a suitable relabelling rule to mma_map
        // otherwise, insert a new rule provided there is no symbol collision
        for(const auto& v : source)
          {
            const GiNaC::ex& label = v.first;
            const GiNaC::ex& value = v.second.subs(subs_rules);
        
            const auto& label_symbol = GiNaC::ex_to<GiNaC::symbol>(label);
        
            // search for an existing remap rule which matches this one
            auto u = std::find_if(dest.begin(), dest.end(), [&](const GiNaC::exmap::value_type& a) -> bool
                                    { return static_cast<bool>(a.second == value); });
        
            // if one exists, just add a relabelling rule that will redirect this Rayleigh symbol
            // to the existing definition
            if(u != dest.end())
              {
                mma_map[label_symbol] = u->first;
                continue;
              }
        
            // otherwise need to add a new relabelling rule
        
            // first, is there a symbol collision?
            bool collision = false;

            if(reserved.find(label_symbol) != reserved.end()) collision = true;
            if(dest.find(label_symbol) != dest.end()) collision = true;
        
            // if no collision, just keep the same symbol to avoid proliferation
            if(!collision)
              {
                dest[label_symbol] = value;
                continue;
              }
        
            // otherwise, need to manufacture a new symbol
            auto relabel = sf.make_unique_Rayleigh_momentum();
        
            mma_map[label_symbol] = relabel;
            dest[relabel] = value;
          }
        
        return mma_map;
      }
    
    
    GiNaC::exmap
    merge_Rayleigh_lists(const GiNaC::exmap& source, GiNaC::exmap& dest, GiNaC_symbol_set& reserved,
                         const GiNaC::exmap& subs_rules, symbol_factory& sf)
      {
        auto Ray_remap = merge_Rayleigh_rules(dest, source, reserved, subs_rules, sf);
    
        // merge symbols in Rayleigh_list into reserved symbols
        std::for_each(dest.begin(), dest.end(), [&](const GiNaC::exmap::value_type& v) -> void
          {
            const GiNaC::ex s = v.first;
            if(GiNaC::is_a<GiNaC::symbol>(s)) reserved.insert(GiNaC::ex_to<GiNaC::symbol>(s));
          });
        
        return Ray_remap;
      }
    
    
    GiNaC::exmap remove_Rayleigh_trivial(GiNaC::exmap& rm)
      {
        GiNaC::exmap trivmap;

        for(auto t = rm.begin(); t != rm.end(); /* deliberately left empty */)
          {
            const GiNaC::ex& rhs = t->second;
            
            // if trivially just a single symbol, perform replacement
            if(GiNaC::is_a<GiNaC::symbol>(rhs))
              {
                trivmap.emplace(t->first, t->second);
                t = rm.erase(t);
                continue;
              }
            
            // alternatively, if just a numerical factor multiplied by a symbol, also perform replacement
            if(GiNaC::is_a<GiNaC::mul>(rhs) && rhs.nops() == 2)
              {
                const auto& rhs_mul = GiNaC::ex_to<GiNaC::mul>(rhs);
                const GiNaC::ex& op1 = rhs_mul.op(0);
                const GiNaC::ex& op2 = rhs_mul.op(1);
                
                if((GiNaC::is_a<GiNaC::symbol>(op1) && GiNaC::is_a<GiNaC::numeric>(op2))
                   || (GiNaC::is_a<GiNaC::symbol>(op2) && GiNaC::is_a<GiNaC::numeric>(op1)))
                  {
                    trivmap.emplace(t->first, t->second);
                    t = rm.erase(t);
                    continue;
                  }
              }
            
            ++t;
          }
        
        return trivmap;
      }
    
  }   // namespace detail
