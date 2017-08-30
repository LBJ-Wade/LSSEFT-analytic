//
// Created by David Seery on 30/08/2017.
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

#include "Rayleigh_momenta.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


namespace detail
  {
    
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
            
//            if(static_cast<bool>(value == 0))
//              throw exception(ERROR_RAYLEIGH_MOMENTUM_IS_ZERO, exception_code::Rayleigh_error);
            
            const auto& label_symbol = GiNaC::ex_to<GiNaC::symbol>(label);
            
            // search for an existing remap rule which matches this one (or its negative)
            auto u = std::find_if(dest.begin(), dest.end(), [&](const GiNaC::exmap::value_type& a) -> bool
              { return static_cast<bool>(a.second == value); });
            
            // if one exists, just add a relabelling rule that will redirect this Rayleigh symbol
            // to the existing definition
            if(u != dest.end())
              {
                mma_map[label_symbol] = u->first;
                continue;
              }

            u = std::find_if(dest.begin(), dest.end(), [&](const GiNaC::exmap::value_type& a) -> bool
              { return static_cast<bool>(a.second == -value); });

            // if one exists, just add a relabelling rule that will redirect this Rayleigh symbol
            // to the existing definition
            if(u != dest.end())
              {
                mma_map[label_symbol] = -u->first;
                continue;
              }
            
            // there was no match, so need to add a new relabelling rule
            
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
    
    
    void prune_Rayleigh_list(GiNaC::exmap& list, const GiNaC::ex& K)
      {
        // get symbols used in K
        const auto& symbols = get_expr_symbols(K);
        
        // step through 'list' removing replacement rules for symbols that no longer occur in K
        for(auto t = list.begin(); t != list.end(); /* intentionally left blank*/)
          {
            const auto& sym = GiNaC::ex_to<GiNaC::symbol>(t->first);
            
            // if symbols occurs in K, move on
            if(symbols.find(sym) != symbols.end())
              {
                ++t;
                continue;
              }
            
            // remove this replacement rule
            t = list.erase(t);
          }
      }
    
  }   // namespace detail
