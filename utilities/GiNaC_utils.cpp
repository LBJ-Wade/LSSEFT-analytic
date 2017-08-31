//
// Created by David Seery on 22/08/2017.
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

#include <localizations/messages.h>
#include <shared/exceptions.h>
#include "GiNaC_utils.h"


bool is_rational(const GiNaC::ex& expr, GiNaC::exvector dummies);


GiNaC_symbol_set get_expr_symbols(const GiNaC::ex& expr)
  {
    GiNaC_symbol_set syms;
    
    if(GiNaC::is_exactly_a<GiNaC::indexed>(expr))
      {
        return get_expr_symbols(expr.op(0));
      }
    
    if(GiNaC::is_exactly_a<GiNaC::symbol>(expr))
      {
        syms.insert(GiNaC::ex_to<GiNaC::symbol>(expr));
        return syms;
      }

    size_t nops = expr.nops();
    for(size_t i = 0; i < nops; ++i)
      {
        auto new_syms = get_expr_symbols(expr.op(i));
        std::copy(new_syms.begin(), new_syms.end(), std::inserter(syms, syms.end()));
      }
    
    return syms;
  }


void get_expr_indices_impl(const GiNaC::ex& expr, std::map<GiNaC::symbol, size_t>& sym_list)
  {
    if(GiNaC::is_exactly_a<GiNaC::indexed>(expr))
      {
        // walk through all indices, extracting their symbols
        size_t nops = expr.nops();
        for(size_t i = 1; i < nops; ++i)
          {
            const auto& idx_raw = expr.op(i);
            if(GiNaC::is_a<GiNaC::idx>(idx_raw))
              {
                const auto& idx = GiNaC::ex_to<GiNaC::idx>(idx_raw);
                const auto& sym = idx.get_value();
                if(GiNaC::is_exactly_a<GiNaC::symbol>(sym))
                  {
                    sym_list[GiNaC::ex_to<GiNaC::symbol>(sym)]++; // safe since values default to zero
                  }
              }
          }
      }
    
    size_t nops = expr.nops();
    for(size_t i = 0; i < nops; ++i)
      {
        get_expr_indices_impl(expr.op(i), sym_list);
      }
  }


GiNaC_symbol_set get_expr_indices(const GiNaC::ex& expr, size_t min_occurrences)
  {
    using symmap = std::map< GiNaC::symbol, size_t >;
    symmap list;
    get_expr_indices_impl(expr, list);

    GiNaC_symbol_set idxs;
    std::for_each(list.begin(), list.end(), [&](const symmap::value_type& v) -> void
      {
        if(v.second >= min_occurrences) idxs.insert(v.first);
      });
    
    return idxs;
  }


GiNaC::ex simplify_add(const GiNaC::ex& expr, const GiNaC::scalar_products& sp)
  {
    GiNaC::ex val{0};
    
    for(auto t = expr.begin(); t != expr.end(); ++t)
      {
        val += simplify_index(*t, sp);
      }
    
    return val;
  }


GiNaC::ex simplify_mul(const GiNaC::ex& expr, const GiNaC::scalar_products& sp)
  {
    GiNaC::ex val{1};
    
    for(auto t = expr.begin(); t != expr.end(); ++t)
      {
        val *= simplify_index(*t, sp);
      }
    
    return val.simplify_indexed(sp);
  }


GiNaC::ex simplify_pow(const GiNaC::ex& expr, const GiNaC::scalar_products& sp)
  {
    const GiNaC::ex& base = expr.op(0);
    const GiNaC::ex& exponent = expr.op(1);
    
    // if the base is not indexed then apply recursively and return
    if(!GiNaC::is_exactly_a<GiNaC::indexed>(base))
      return GiNaC::pow(simplify_index(base.expand(GiNaC::expand_options::expand_indexed), sp), exponent);

    // if base has too many indices then give up
    if(base.nops() > 2) return expr;

    // do nothing if the exponent isn't a number; the expression probably doesn't make sense
    if(!GiNaC::is_exactly_a<GiNaC::numeric>(exponent)) return expr;
    
    const auto& exp_as_numeric = GiNaC::ex_to<GiNaC::numeric>(exponent);
    
    // do nothing if the exponent isn't an integer; the expression probably doesn't make sense
    if(!GiNaC::is_integer(exp_as_numeric)) return expr;
    
    const int p = exp_as_numeric.to_int();
    
    // do nothing if the exponent isn't an even number; the expression probably doesn't make sense
    if(std::abs(p) % 2 != 0) return expr;

    return GiNaC::pow((base*base).simplify_indexed(sp), p/2);
  }


GiNaC::ex simplify_index(const GiNaC::ex& expr, const GiNaC::scalar_products& sp)
  {
    if(GiNaC::is_exactly_a<GiNaC::add>(expr))
      {
        return simplify_add(expr, sp);
      }
    if(GiNaC::is_exactly_a<GiNaC::mul>(expr))
      {
        return simplify_mul(expr, sp);
      }
    if(GiNaC::is_exactly_a<GiNaC::power>(expr))
      {
        return simplify_pow(expr, sp);
      }
    if(GiNaC::is_exactly_a<GiNaC::indexed>(expr))
      {
        // exactly an indexed object,
        return expr.simplify_indexed(sp);
      }
    
    // nothing to do, so return expression unaltered
    return expr;
  }


GiNaC::ex simplify_index(const GiNaC::ex& expr)
  {
    return simplify_index(expr.expand(GiNaC::expand_options::expand_indexed), GiNaC::scalar_products{});
  }


bool is_rational_add(const GiNaC::ex& expr)
  {
    for(auto t = expr.begin(); t != expr.end(); ++t)
      {
        if(!is_rational(*t)) return false;
      }
    
    return true;
  }


bool is_rational_mul(const GiNaC::ex& expr)
  {
    GiNaC::exvector dummies = GiNaC::get_all_dummy_indices_safely(expr);

    for(auto t = expr.begin(); t != expr.end(); ++t)
      {
        if(!is_rational(*t, dummies)) return false;
      }

    return true;
  }


bool is_rational_pow(const GiNaC::ex& expr)
  {
    const GiNaC::ex& base = expr.op(0);
    const GiNaC::ex& exponent = expr.op(1);
    
    // if the power is not a number then this is not a rational function
    if(!GiNaC::is_exactly_a<GiNaC::numeric>(exponent)) return false;
    
    const auto& exp_as_numeric = GiNaC::ex_to<GiNaC::numeric>(exponent);
    
    // if the power is not an integer then this is not a rational function
    if(!GiNaC::is_integer(exp_as_numeric)) return false;
    
    // if the base is not an indexed expression, then the power is rational if the base is
    if(!GiNaC::is_exactly_a<GiNaC::indexed>(base))
      return is_rational(base.expand(GiNaC::expand_options::expand_indexed));
    
    // otherwise, this is rational if the power is divisible by two
    int p = exp_as_numeric.to_int();
    return (std::abs(p) % 2) == 0;
  }


bool is_rational(const GiNaC::ex& expr, GiNaC::exvector dummies)
  {
    if(GiNaC::is_exactly_a<GiNaC::add>(expr))
      {
        return is_rational_add(expr);
      }
    if(GiNaC::is_exactly_a<GiNaC::mul>(expr))
      {
        return is_rational_mul(expr);
      }
    if(GiNaC::is_exactly_a<GiNaC::power>(expr))
      {
        return is_rational_pow(expr);
      }
    if(GiNaC::is_exactly_a<GiNaC::indexed>(expr))
      {
        // get indices on this object
        GiNaC::exvector indices = expr.get_free_indices();
    
        // determine which of these are not dummy indices
        GiNaC::exvector not_dummies;
        std::set_difference(indices.begin(), indices.end(), dummies.begin(), dummies.end(),
                            std::back_inserter(not_dummies));
        
        // object is not rational if there is a non-dummy index
        return not_dummies.empty();
      }

    return GiNaC::is_exactly_a<GiNaC::symbol>(expr)
           || GiNaC::is_exactly_a<GiNaC::numeric>(expr);
  }


bool is_rational(const GiNaC::ex& expr)
  {
    return is_rational(expr.expand(GiNaC::expand_options::expand_indexed), GiNaC::exvector{});
  }


GiNaC::exvector to_exvector(const GiNaC::ex& expr)
  {
    GiNaC::exvector vec;
    vec.reserve(expr.nops());

    for(size_t i = 0; i < expr.nops(); ++i)
      {
        vec.push_back(expr.op(i));
      }

    return vec;
  }
