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

#include "loop_integral.h"

#include "detail/special_functions.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


std::ostream& operator<<(std::ostream& str, const loop_integral& obj)
  {
    obj.write(str);
    return str;
  }


void loop_integral::write(std::ostream& out) const
  {
    std::cout << "  time function = " << this->tm << '\n';
    std::cout << "  momentum kernel = " << this->K << '\n';
    std::cout << "  Wick product = " << this->WickString << '\n';

    std::cout << "  loop momenta =";
    for(const auto& sym : this->loop_momenta)
      {
        std::cout << " " << sym;
      }
    if(this->loop_momenta.empty()) std::cout << " <none>";
    std::cout << '\n';

    std::cout << "  Rayleigh momenta =";
    for(const auto& rule : this->Rayleigh_momenta)
      {
        std::cout << " " << rule.first << " -> " << rule.second << ";";
      }
    if(this->Rayleigh_momenta.empty()) std::cout << " <none>";
    std::cout << '\n';
  }


namespace inner_products_to_cos_impl
  {

    GiNaC::ex convert(const GiNaC::ex& expr);


    GiNaC::ex convert_add(const GiNaC::ex& expr)
      {
        GiNaC::ex val{0};

        for(auto t = expr.begin(); t != expr.end(); ++t)
          {
            val += convert(*t);
          }

        return val;
      }


    GiNaC::ex convert_mul(const GiNaC::ex& expr)
      {
        auto exvec = to_exvector(expr);
        GiNaC::ex val{1};

        for(auto t = exvec.begin(); t != exvec.end(); /* intentionally left blank*/)
          {
            // does this factor carry an index? if not push it to the result and carry on
            if(!GiNaC::is_a<GiNaC::indexed>(*t))
              {
                val *= convert(*t);
                ++t;
                continue;
              }

            const auto& idx_item = GiNaC::ex_to<GiNaC::indexed>(*t);
            if(idx_item.nops() > 2)
              throw exception(ERROR_CANT_HANDLE_TENSORS, exception_code::loop_transformation_error);

            const auto idx_cand = idx_item.op(1);
            if(!GiNaC::is_a<GiNaC::idx>(idx_cand))
              throw exception(ERROR_EXPECTED_INDEX_LABEL, exception_code::loop_transformation_error);

            const auto& idx_label = GiNaC::ex_to<GiNaC::idx>(idx_cand);

            // now find the partner object among the other factors
            auto u = t+1;
            for(; u != exvec.end(); ++u)
              {
                if(!GiNaC::is_a<GiNaC::indexed>(*u)) continue;

                const auto idx_item2 = GiNaC::ex_to<GiNaC::indexed>(*u);
                if(idx_item2.nops() > 2)
                  throw exception(ERROR_CANT_HANDLE_TENSORS, exception_code::loop_transformation_error);

                const auto idx_cand2 = idx_item2.op(1);
                if(!GiNaC::is_a<GiNaC::idx>(idx_cand2))
                  throw exception(ERROR_EXPECTED_INDEX_LABEL, exception_code::loop_transformation_error);

                const auto& idx_label2 = GiNaC::ex_to<GiNaC::idx>(idx_cand2);

                if(idx_label2 == idx_label) break;
              }

            if(u == exvec.end())
              throw exception(ERROR_COULD_NOT_FIND_PARTNER_INDEX, exception_code::loop_transformation_error);

            // dot product is between base elements
            const auto base1 = idx_item.op(0);
            if(!GiNaC::is_a<GiNaC::symbol>(base1))
              throw exception(ERROR_EXPECTED_SYMBOL, exception_code::loop_transformation_error);

            auto sym1 = GiNaC::ex_to<GiNaC::symbol>(base1);

            const auto base2 = u->op(0);
            if(!GiNaC::is_a<GiNaC::symbol>(base2))
              throw exception(ERROR_EXPECTED_SYMBOL, exception_code::loop_transformation_error);

            auto sym2 = GiNaC::ex_to<GiNaC::symbol>(base2);

            if(std::less<GiNaC::symbol>{}(sym2, sym1)) std::swap(sym1, sym2);

            // replace with Cos factor, if needed
            if(static_cast<bool>(base1 == base2))
              {
                val *= sym1 * sym1;
              }
            else
              {
                val *= sym1 * sym2 * Angular::Cos(sym1, sym2);
              }

            // erase factors that have been replaced
            exvec.erase(u);   // only invalidates iterators from u onwards
            t = exvec.erase(t);
          }

        return val;
      }


    GiNaC::ex convert_power(const GiNaC::ex expr)
      {
        const GiNaC::ex& base = expr.op(0);
        const GiNaC::ex& exponent = expr.op(1);

        // if the base is not an indexed object then apply recursively
        if(!GiNaC::is_a<GiNaC::indexed>(base))
          return GiNaC::pow(convert(base), exponent);

        // if base has too many indices then complain
        if(base.nops() > 2)
          throw exception(ERROR_CANT_HANDLE_TENSORS, exception_code::loop_transformation_error);

        // get exponent as a numeric
        if(!GiNaC::is_a<GiNaC::numeric>(exponent))
          throw exception(EXPECTED_EXPONENT_TO_BE_NUMERIC, exception_code::loop_transformation_error);

        const auto& exp_as_num = GiNaC::ex_to<GiNaC::numeric>(exponent);

        // if exponent isn't an integer then complain
        if(!GiNaC::is_integer(exp_as_num))
          throw exception(EXPECTED_EXPONENT_TO_BE_INTEGER, exception_code::loop_transformation_error);

        const int p = exp_as_num.to_int();

        if(std::abs(p) % 2 != 0) throw exception(EXPECTED_EXPONENT_TO_BE_EVEN_INTEGER, exception_code::loop_transformation_error);

        // strip out kernel from base
        const GiNaC::ex& symb = base.op(0);

        if(!GiNaC::is_a<GiNaC::symbol>(symb))
          throw exception(ERROR_EXPECTED_SYMBOL, exception_code::loop_transformation_error);

        return GiNaC::pow(symb*symb, p/2);
      }


    GiNaC::ex convert(const GiNaC::ex& expr)
      {
        auto expr_expand = expr.expand(GiNaC::expand_options::expand_indexed);

        if(GiNaC::is_a<GiNaC::add>(expr_expand))
          {
            return convert_add(expr_expand);
          }
        if(GiNaC::is_a<GiNaC::mul>(expr_expand))
          {
            return convert_mul(expr_expand);
          }
        if(GiNaC::is_a<GiNaC::power>(expr_expand))
          {
            return convert_power(expr_expand);
          }

        // otherwise assume nothing to do, so return
        return expr_expand;
      }

  }   // namespace inner_products_to_cos_impl


void loop_integral::inner_products_to_cos()
  {
    // step through the kernel, converting any explicit inner products into cosines;
    // these can later be converted into Legendre polynomials if required
    GiNaC::ex new_K = inner_products_to_cos_impl::convert(this->K);
    this->K = new_K;
  }


void loop_integral::canonicalize_momenta()
  {
    GiNaC_symbol_set new_momenta;
    GiNaC::exmap relabel;

    unsigned int count = 0;
    for(const auto& l : this->loop_momenta)
      {
        const auto L = this->sf.make_canonical_loop_momentum(count++);
        relabel[l] = L;
        new_momenta.insert(L);
      }

    this->K = this->K.subs(relabel);
    this->WickString = this->WickString.subs(relabel);

    for(auto& rule : this->Rayleigh_momenta)
      {
        rule.second = rule.second.subs(relabel);
      }

    this->loop_momenta = new_momenta;
  }
