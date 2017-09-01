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
      : integrand(ig_),
        measure(ms_),
        WickProduct(wp_),
        tm(tm_),
        variables(vs_)
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
  : loop_int(i_)
  {
    if(loop_int.get_loop_order() == 0)
      throw exception(ERROR_ONE_LOOP_REDUCE_WITH_TREE, exception_code::loop_transformation_error);

    if(loop_int.get_loop_order() > 1)
      throw exception(ERROR_ONE_LOOP_REDUCE_WITH_MULTIPLE_LOOPS, exception_code::loop_transformation_error);

    // extract momentum kernel from integral
    auto K = loop_int.get_kernel();
    K = dot_products_to_cos(K);
    K = K.expand();

    // apply term-by-term decomposition to K
    if(GiNaC::is_a<GiNaC::mul>(K))
      {
        // just a single product at the top level
        this->reduce(K);
      }
    else if(GiNaC::is_a<GiNaC::add>(K))
      {
        // a sum of terms at top level; work through each one
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
    // get Rayleigh momenta from integral
    const auto& Rayleigh_momenta = this->loop_int.get_Rayleigh_momenta();

    // find which Rayleigh momenta this term depends on, if any
    GiNaC_symbol_set Rayleigh_mma;

    for(const auto& rule : Rayleigh_momenta)
      {
        const auto& sym = GiNaC::ex_to<GiNaC::symbol>(rule.first);
        if(term.has(sym))
          {
            Rayleigh_mma.insert(sym);
          }
      }

    // if no Rayleigh momenta then we need only account for the loop momentum
    if(Rayleigh_mma.empty())
      {
        this->one_loop_reduce_zero_Rayleigh(term);
        return;
      }

    // if exactly one Rayleigh momentum then we know how to reduce it
    if(Rayleigh_mma.size() == 1)
      {
        this->one_loop_reduce_one_Rayleigh(term, *Rayleigh_mma.begin());
        return;
      }

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
    temp = cosines_to_Legendre(term, L);

    // construct integration element
    auto measure = L*L / GiNaC::pow(2*GiNaC::Pi, 3);

    using one_loop_reduced_integral_impl::integration_element;
    using one_loop_reduced_integral_impl::key;
    auto elt = std::make_unique<integration_element>(temp, measure, WickProduct, tm, GiNaC_symbol_set{L});

    // insert in database
    this->integrand[key{*elt}].push_back(std::move(elt));
  }


void one_loop_reduced_integral::one_loop_reduce_one_Rayleigh(const GiNaC::ex& term, const GiNaC::symbol& R)
  {
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
