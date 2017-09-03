//
// Created by David Seery on 01/09/2017.
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

#include "one_loop_reduced_integral.h"

#include "detail/special_functions.h"
#include "detail/legendre_utils.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


namespace one_loop_reduced_integral_impl
  {

    integration_element::integration_element(GiNaC::ex ig_, GiNaC::ex ms_, GiNaC::ex wp_, time_function tm_,
                                             GiNaC_symbol_set vs_)
      : integrand(std::move(ig_)),
        measure(std::move(ms_)),
        WickProduct(std::move(wp_)),
        tm(std::move(tm_)),
        variables(std::move(vs_))
      {
      }


    void integration_element::write(std::ostream& str) const
      {
        str << "integral";
        for(const auto& sym : this->variables)
          {
            str << " d" << sym;
          }
        str << '\n';

        str << "  time function = " << this->tm << '\n';
        str << "  measure = " << this->measure << '\n';
        str << "  Wick product = " << this->WickProduct << '\n';
        str << "  integrand = " << this->integrand << '\n';
      }


    key::key(const integration_element& elt_)
      : elt(elt_)
      {
      }


    size_t key::hash() const
      {
        // concatenate symbols in integration variable list
        std::string symbol_string;
        std::for_each(this->elt.variables.begin(), this->elt.variables.end(),
                      [&](const GiNaC::symbol& e) -> std::string
                        { return symbol_string += e.get_name(); });

        // combine both hashes together
        size_t h = 0;
        hash_impl::hash_combine(h, symbol_string);

        // return final value
        return h;
      }


    bool key::is_equal(const key& obj) const
      {
        return std::equal(this->elt.variables.cbegin(), this->elt.variables.cend(),
                          obj.elt.variables.cbegin(), obj.elt.variables.cend(),
                          [](const GiNaC::symbol& asym, const GiNaC::symbol& bsym) -> bool
                            { return asym.get_name() == bsym.get_name(); });
      }


    std::ostream& operator<<(std::ostream& str, const integration_element& obj)
      {
        obj.write(str);
        return str;
      }

  }   // namespace one_loop_reduced_integral_impl


one_loop_reduced_integral::one_loop_reduced_integral(const loop_integral& i_)
  : loop_int(i_),
    Rayleigh_momenta(i_.get_Rayleigh_momenta())
  {
    // throw if we were given a tree-level expression
    if(loop_int.get_loop_order() == 0)
      throw exception(ERROR_ONE_LOOP_REDUCE_WITH_TREE, exception_code::loop_transformation_error);

    // throw if we were given a 2+ loop expression
    if(loop_int.get_loop_order() > 1)
      throw exception(ERROR_ONE_LOOP_REDUCE_WITH_MULTIPLE_LOOPS, exception_code::loop_transformation_error);

    // extract momentum kernel from integral
    auto K = loop_int.get_kernel();

    // convert explicit dot products to Cos(a,b) format and expand to get a representation
    // suitable for term-by-term decomposition into a sum of products of Legendre polynomials
    K = dot_products_to_cos(K);
    K = K.expand();

    // apply term-by-term decomposition to K
    if(GiNaC::is_a<GiNaC::mul>(K))
      {
        // just a single product at the top level, so reduce that
        this->reduce(K);
      }
    else if(GiNaC::is_a<GiNaC::add>(K))
      {
        // a sum of terms at top level; work through each one and reduce them term-by-term
        for(size_t i = 0; i < K.nops(); ++i)
          {
            this->reduce(K.op(i));
          }
      }
    else
      throw exception(ERROR_BADLY_FORMED_TOP_LEVEL_MOMENTUM_KERNEL, exception_code::loop_transformation_error);
  }


void one_loop_reduced_integral::reduce(const GiNaC::ex& term)
  {
    // find which Rayleigh momenta this term depends on, if any
    GiNaC_symbol_set Rayleigh_mma;

    for(const auto& rule : this->Rayleigh_momenta)
      {
        const auto& sym = GiNaC::ex_to<GiNaC::symbol>(rule.first);
        if(term.has(sym))
          {
            Rayleigh_mma.insert(sym);
          }
      }

    if(Rayleigh_mma.empty())
      {
        // if no Rayleigh momenta then we need only account for the loop momentum
        this->one_loop_reduce_zero_Rayleigh(term);
      }
    else if(Rayleigh_mma.size() == 1)
      {
        // if exactly one Rayleigh momentum then we know how to reduce it
        this->one_loop_reduce_one_Rayleigh(term, *Rayleigh_mma.begin());
      }
    else
      // currently we don't know how to reduce terms with more than one Rayleigh momentum
      throw exception(ERROR_MULTIPLE_RAYLEIGH_MOMENTA_NOT_IMPLEMENTED, exception_code::loop_transformation_error);
  }


void one_loop_reduced_integral::one_loop_reduce_zero_Rayleigh(const GiNaC::ex& term)
  {
    // extract loop momentum from integral
    const auto& lm = this->loop_int.get_loop_momenta();
    const auto& L = *lm.begin();

    // extract Wick product from integral
    const auto& WickProduct = this->loop_int.get_Wick_product();

    // extract time factor from integral
    const auto& tm = this->loop_int.get_time_function();

    // first, convert all angular terms involving the loop momentum to Legendre representation
    auto temp = Legendre_to_cosines(term, L);
    temp = cosines_to_Legendre(temp, L);

    // step through expression, identifying terms with zero, one, two or more Legendre polynomials
    temp = this->integrate_Legendre(temp, L);

    if(temp != 0)
      {
        // construct integration element
        auto measure = L*L / GiNaC::pow(2*GiNaC::Pi, 3);

        using one_loop_reduced_integral_impl::integration_element;
        using one_loop_reduced_integral_impl::key;
        auto elt = std::make_unique<integration_element>(temp, measure, WickProduct, tm, GiNaC_symbol_set{L});

        // insert in database
        this->integrand[key{*elt}].push_back(std::move(elt));
      }
  }


void one_loop_reduced_integral::one_loop_reduce_one_Rayleigh(const GiNaC::ex& term, const GiNaC::symbol& R)
  {
  }


GiNaC::ex one_loop_reduced_integral::integrate_Legendre(const GiNaC::ex& term, const GiNaC::ex& q)
  {
    if(GiNaC::is_a<GiNaC::mul>(term))
      {
        return this->apply_Legendre_orthogonality(term, q);
      }

    if(GiNaC::is_a<GiNaC::add>(term))
      {
        GiNaC::ex temp{0};
        for(size_t i = 0; i < temp.nops(); ++i)
          {
            temp += this->apply_Legendre_orthogonality(term.op(i), q);
          }
        return temp;
      }

    throw exception(ERROR_BADLY_FORMED_TOP_LEVEL_LEGENDRE_SUM, exception_code::loop_transformation_error);
  }


GiNaC::ex one_loop_reduced_integral::apply_Legendre_orthogonality(const GiNaC::ex& expr, const GiNaC::ex& q)
  {
    GiNaC::ex temp{1};

    if(!GiNaC::is_a<GiNaC::mul>(expr))
      throw exception(ERROR_BADLY_FORED_LEGENDRE_SUM_TERM, exception_code::loop_transformation_error);

    using Legendre_list = std::vector< std::pair< GiNaC::symbol, unsigned int > >;
    Legendre_list partner_q;

    for(size_t i = 0; i < expr.nops(); ++i)
      {
        const auto term = expr.op(i);
        if(!GiNaC::is_a<GiNaC::function>(term)) { temp *= term; continue; }

        const auto& fn = GiNaC::ex_to<GiNaC::function>(term);
        if(fn.get_name() != "LegP") { temp *= term; continue; }

        auto n = static_cast<unsigned int>(GiNaC::ex_to<GiNaC::numeric>(fn.op(0)).to_int());
        auto p1 = GiNaC::ex_to<GiNaC::symbol>(fn.op(1));
        auto p2 = GiNaC::ex_to<GiNaC::symbol>(fn.op(2));

        if(p1 == q) { partner_q.emplace_back(p2, n); continue; }
        if(p2 == q) { partner_q.emplace_back(p1, n); continue; }

        // not a Legendre polynomial involving q, so move on
        temp *= term;
      }

    // if no Legendre polynomials, equivalent to LegP(0, x)
    if(partner_q.empty())
      {
        return GiNaC::numeric{4} * GiNaC::Pi * temp;
      }

    // if one Legendre polynomial, gives nonzero if n=0
    if(partner_q.size() == 1)
      {
        if(partner_q.front().second != 0) return GiNaC::ex{0};
        return GiNaC::numeric{4} * GiNaC::Pi * temp;
      }

    // if more than two Legendre polynomials, don't know what to do
    if(partner_q.size() > 2)
        throw exception(ERROR_CANT_INTEGRATE_MORE_THAN_TWO_LEGP, exception_code::loop_transformation_error);

    if(partner_q.front().second != partner_q.back().second) return GiNaC::ex{0};

    unsigned int n = partner_q.front().second;
    auto p1 = partner_q.front().first;
    auto p2 = partner_q.back().first;

    if(std::less<GiNaC::symbol>{}(p2, p1)) std::swap(p1, p2);

    return GiNaC::numeric{4} * GiNaC::Pi / (2*GiNaC::numeric{n} + 1) * Angular::LegP(n, p1, p2) * temp;
  }


void one_loop_reduced_integral::write(std::ostream& out) const
  {
    using one_loop_reduced_integral_impl::integration_element;

    for(const auto& record : this->integrand)
      {
        const auto& vars = record.first;
        const auto& data = record.second;

        for(const auto& iptr : data)
          {
            const auto& i = *iptr;
            out << i;
          }
      }
  }


std::ostream& operator<<(std::ostream& str, const one_loop_reduced_integral& obj)
  {
    obj.write(str);
    return str;
  }
