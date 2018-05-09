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

#include "GiNaC_utils.h"

#include "services/service_locator.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


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
            const auto idx_raw = expr.op(i);
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


// forward declare main simplify_index implementation method
GiNaC::ex simplify_index_impl(const GiNaC::ex& expr, const GiNaC::scalar_products& sp, const GiNaC::exmap& Rayleigh_list,
                              service_locator& loc);


GiNaC::ex simplify_add(const GiNaC::ex& expr, const GiNaC::scalar_products& sp, const GiNaC::exmap& Rayleigh_list,
                       service_locator& loc)
  {
    GiNaC::ex val{0};
    
    for(auto t = expr.begin(); t != expr.end(); ++t)
      {
        val += simplify_index_impl(*t, sp, Rayleigh_list, loc);
      }

    // call simplify_indexed() to try to perform index relabelling, if that helps reduce the number of terms
    return val;
  }


std::pair<bool, GiNaC::ex> try_match(const GiNaC::ex& expr, const GiNaC::ex& LHS_search, GiNaC::exmap& LHS_data,
const GiNaC::ex& RHS_search, GiNaC::exmap& RHS_data)
  {
    bool found_LHS = false;
    bool found_RHS = false;

    const auto& term = GiNaC::ex_to<GiNaC::mul>(expr);

    GiNaC::ex val{1};

    for(size_t i = 0; i < term.nops(); ++i)
      {
        const auto& factor = term.op(i);
        bool matched = false;

        if(!found_LHS)
          {
            matched = found_LHS = factor.match(LHS_search, LHS_data);
          }
        if(!matched && !found_RHS)
          {
            matched = found_RHS = factor.match(RHS_search, RHS_data);
          }

        if(!matched) val *= factor;
      }

    return std::make_pair(found_LHS && found_RHS, val);
  };


GiNaC::ex simplify_mul(const GiNaC::ex& expr, const GiNaC::scalar_products& sp, const GiNaC::exmap& Rayleigh_list,
                       service_locator& loc)
  {
    GiNaC::ex val{1};
    
    for(auto t = expr.begin(); t != expr.end(); ++t)
      {
        auto factor = simplify_index_impl(*t, sp, Rayleigh_list, loc);

        val *= factor;
      }

    // try to simplify powers of Rayleigh momenta that also occur in the numerator
    for(const auto& rule : Rayleigh_list)
      {
        bool flag = true;

        while(flag)
          {
            const auto& LHS = rule.first;
            const auto& RHS = rule.second;

            // see if a matching powers of a Rayleigh momentum and its value are present
            GiNaC::exmap LHS_data;
            GiNaC::ex LHS_search = GiNaC::pow(GiNaC::indexed(LHS, loc.get_symbol_factory().make_wildcard_index(0)), GiNaC::wild(1));

            GiNaC::exmap RHS_data;
            GiNaC::ex RHS_search = GiNaC::pow(GiNaC::indexed(RHS, loc.get_symbol_factory().make_wildcard_index(0)), GiNaC::wild(1));

            GiNaC::ex residual;
            std::tie(flag, residual) = try_match(val, LHS_search, LHS_data, RHS_search, RHS_data);

            if(flag)
              {
                const auto& LHS_power = GiNaC::ex_to<GiNaC::numeric>(LHS_data[GiNaC::wild(1)]);
                const auto& RHS_power = GiNaC::ex_to<GiNaC::numeric>(RHS_data[GiNaC::wild(1)]);

                if(LHS_power > 0)
                  throw exception(ERROR_RAYLEIGH_MOMENTA_POSITIVE_POWER, exception_code::Rayleigh_error);

                if(LHS_power.to_int() + RHS_power.to_int() < 0) // more powers in the denominator than the numerator
                  {
                    const auto& index = LHS_data[GiNaC::wild(0)];
                    val = residual * GiNaC::pow(GiNaC::indexed(LHS, index), LHS_power + RHS_power);
                  }
                else if(LHS_power.to_int() + RHS_power.to_int() > 0) // more powers in numerator than denomiantor
                  {
                    const auto& index = RHS_data[GiNaC::wild(0)];
                    val = residual * GiNaC::pow(GiNaC::indexed(RHS, index), LHS_power + RHS_power);
                  }
                else
                  {
                    val = residual;
                  }
              }
          }
      }

    return val;
  }


GiNaC::ex simplify_pow_impl(const GiNaC::ex& base_kernel, const GiNaC::idx& base_idx, const GiNaC::ex& exponent,
                            const GiNaC::scalar_products& sp, const GiNaC::exmap& Rayleigh_list, service_locator& loc)
  {
    auto base = GiNaC::indexed(base_kernel, base_idx);

    // throw if the exponent isn't a number; the expression probably doesn't make sense
    if(!GiNaC::is_exactly_a<GiNaC::numeric>(exponent))
      throw exception(ERROR_EXPONENT_OF_INDEXED_POWER_IS_NOT_NUMERIC, exception_code::index_error);

    const auto& exp_as_numeric = GiNaC::ex_to<GiNaC::numeric>(exponent);

    // throw if the exponent isn't an integer; the expression probably doesn't make sense
    if(!GiNaC::is_integer(exp_as_numeric))
      throw exception(ERROR_EXPONENT_OF_INDEXED_POWER_IS_NOT_INTEGER, exception_code::index_error);

    const int p = exp_as_numeric.to_int();

    // throw if the exponent isn't an even number; the expression probably doesn't make sense
    if(std::abs(p) % 2 != 0)
      throw exception(ERROR_EXPONENT_OF_INDEXED_POWER_IS_NOT_EVEN, exception_code::index_error);

    // if the power is positive, try to replace with a scalar product and return
    // if we can't replace by a scalar then no harm is done
    if(p >= 0) return GiNaC::pow((base*base).simplify_indexed(sp), p/2);

    // if p is negative then we should replace this indexed object by a scalar product if possible, and if not
    // we should write in using a Rayleigh momentum

    // try to replace power with a scalar product
    auto new_base = (base*base).simplify_indexed(sp);

    // check whether new_base has any remaining occurrences of indices; if not, return
    auto idxs = get_expr_indices(new_base);
    if(idxs.empty()) return GiNaC::pow(new_base, p/2);

    // try to find base as the RHS of a rule in Rayleigh list
    auto it = std::find_if(Rayleigh_list.begin(), Rayleigh_list.end(),
                           [&](const GiNaC::exmap::value_type& rule) -> bool
                             {
                               const GiNaC::ex& RHS = rule.second;
                               return static_cast<bool>(RHS == base_kernel);
                             });

    if(it != Rayleigh_list.end()) return GiNaC::pow(GiNaC::indexed(it->first, base_idx), p);

    it = std::find_if(Rayleigh_list.begin(), Rayleigh_list.end(),
                      [&](const GiNaC::exmap::value_type& rule) -> bool
                        {
                          const GiNaC::ex& RHS = rule.second;
                          return static_cast<bool>(RHS == -base_kernel);
                        });

    if(it != Rayleigh_list.end()) return GiNaC::pow(-GiNaC::indexed(it->first, base_idx), p);

    // assume is ok even though we failed to match a Rayleigh momentum
    // (eg. sometimes base is a loop momentum)
    return GiNaC::pow(base, p);
  }


GiNaC::ex simplify_pow(const GiNaC::ex& expr, const GiNaC::scalar_products& sp, const GiNaC::exmap& Rayleigh_list,
                       service_locator& loc)
  {
    const GiNaC::ex& base = expr.op(0);
    const GiNaC::ex& exponent = expr.op(1);

    // if the base is not indexed then apply recursively and return

    if(GiNaC::is_a<GiNaC::indexed>(base))
      {
        // if base has too many indices then give up
        if(base.nops() > 2)
          throw exception(ERROR_CANT_HANDLE_TENSORS, exception_code::index_error);

        const GiNaC::ex& base_kernel = base.op(0);
        const auto& base_idx = GiNaC::ex_to<GiNaC::idx>(base.op(1));

        return simplify_pow_impl(base_kernel, base_idx, exponent, sp, Rayleigh_list, loc);
      }

    if(GiNaC::is_a<GiNaC::add>(base) && base.nops() > 0)
      {
        GiNaC::ex f1 = 1;
        GiNaC::ex t1 = base.op(0);

        // next term in sum could be numerical factor * indexed kernel
        if(GiNaC::is_a<GiNaC::mul>(t1) && t1.nops() == 2)
          {
            if(GiNaC::is_a<GiNaC::numeric>(t1.op(0)) && GiNaC::is_a<GiNaC::indexed>(t1.op(1)))
              {
                f1 = t1.op(0);
                t1 = t1.op(1);
              }
            else if(GiNaC::is_a<GiNaC::numeric>(t1.op(1)) && GiNaC::is_a<GiNaC::indexed>(t1.op(0)))
              {
                f1 = t1.op(1);
                t1 = t1.op(0);
              }
          }

        if(GiNaC::is_a<GiNaC::indexed>(t1))
          {
            if(t1.nops() > 2)
              throw exception(ERROR_CANT_HANDLE_TENSORS, exception_code::index_error);

            const auto& base_idx = GiNaC::ex_to<GiNaC::idx>(t1.op(1));
            auto base_kernel = f1*t1.op(0);

            for(size_t i = 1; i < base.nops(); ++i)
              {
                GiNaC::ex f = 1;
                GiNaC::ex t = base.op(i);

                // next term in sum could be numerical factor * indexed kernel
                if(GiNaC::is_a<GiNaC::mul>(t) && t.nops() == 2)
                  {
                    if(GiNaC::is_a<GiNaC::numeric>(t.op(0)) && GiNaC::is_a<GiNaC::indexed>(t.op(1)))
                      {
                        f = t.op(0);
                        t = t.op(1);
                      }
                    else if(GiNaC::is_a<GiNaC::numeric>(t.op(1)) && GiNaC::is_a<GiNaC::indexed>(t.op(0)))
                      {
                        f = t.op(1);
                        t = t.op(0);
                      }
                  }

                if(!GiNaC::is_a<GiNaC::indexed>(t))
                  {
                    std::cout << base << '\n';
                    throw exception(ERROR_INCONSISTENT_INDEXED_SUM, exception_code::index_error);
                  }

                if(t.nops() > 2)
                  throw exception(ERROR_CANT_HANDLE_TENSORS, exception_code::index_error);

                const GiNaC::ex& t_base = t.op(0);
                const auto& t_idx = GiNaC::ex_to<GiNaC::idx>(t.op(1));

                if(!static_cast<bool>(base_idx == t_idx))
                  {
                    std::cout << base << '\n';
                    std::cout << "base idx = " << base_idx << ", term index = " << t_idx << ", comparison = " << static_cast<bool>(base_idx == t_idx) << '\n';
                    throw exception(ERROR_INCONSISTENT_INDEXED_SUM, exception_code::index_error);
                  }

                base_kernel += f*t_base;
              }

            return simplify_pow_impl(base_kernel, base_idx, exponent, sp, Rayleigh_list, loc);
          }
      }

    return GiNaC::pow(
      simplify_index_impl(
        base.expand(GiNaC::expand_options::expand_indexed), sp, Rayleigh_list, loc
      ), exponent
    );
  }


GiNaC::ex simplify_index_impl(const GiNaC::ex& expr, const GiNaC::scalar_products& sp, const GiNaC::exmap& Rayleigh_list,
                              service_locator& loc)
  {
    if(GiNaC::is_exactly_a<GiNaC::add>(expr))
      {
        return simplify_add(expr, sp, Rayleigh_list, loc);
      }
    if(GiNaC::is_exactly_a<GiNaC::mul>(expr))
      {
        return simplify_mul(expr, sp, Rayleigh_list, loc);
      }
    if(GiNaC::is_exactly_a<GiNaC::power>(expr))
      {
        return simplify_pow(expr, sp, Rayleigh_list, loc);
      }

    // nothing to do, so return expression unaltered
    return expr;
  }


GiNaC::ex simplify_index(const GiNaC::ex& expr, const GiNaC::scalar_products& sp, const GiNaC::exmap& Rayleigh_list,
                         service_locator& loc)
  {
    // filter out any zero elements from Rayleigh list
    GiNaC::exmap new_Rayleigh_list;
    for(const auto& rule : Rayleigh_list)
      {
        const auto& LHS = rule.first;
        const auto& RHS = rule.second;

        if(!static_cast<bool>(RHS == 0)) new_Rayleigh_list[LHS] = RHS;
        else                             new_Rayleigh_list[LHS] = LHS;  // insert identity map
      }

    auto expr_mod = GiNaC::simplify_indexed(expr.expand(GiNaC::expand_options::expand_indexed), sp);

    auto result = simplify_index_impl(expr_mod, sp, new_Rayleigh_list, loc);

    return GiNaC::simplify_indexed(result.expand(GiNaC::expand_options::expand_indexed), sp);
  }


GiNaC::ex simplify_index(const GiNaC::ex& expr, const GiNaC::exmap& Rayleigh_list, service_locator& loc)
  {
    // pass to simplify_index with empty scalar product set
    return simplify_index(expr, GiNaC::scalar_products{}, Rayleigh_list, loc);
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
    const GiNaC::ex base = expr.op(0);
    const GiNaC::ex exponent = expr.op(1);
    
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


std::vector<GiNaC::symbol> order_symbol_set(const GiNaC_symbol_set& syms)
  {
    std::vector<GiNaC::symbol> ordered_set;

    for(const auto& sym : syms)
      {
        ordered_set.push_back(sym);
      }

    std::sort(ordered_set.begin(), ordered_set.end(), std::less<GiNaC::symbol>{});

    return ordered_set;
  }
