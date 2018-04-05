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

#include <string>
#include <sstream>

#include "tree_expression.h"

#include "utilities/hash_combine.h"
#include "utilities/GiNaC_print.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


std::ostream& operator<<(std::ostream& str, const tree_expression& obj)
  {
    obj.write(str);
    return str;
  }


tree_expression::tree_expression(time_function tm_, GiNaC::ex K_, GiNaC::ex ws_, GiNaC_symbol_set em_,
                                 service_locator& lc_)
  : tm(std::move(tm_)),
    K(std::move(K_)),
    WickProduct(std::move(ws_)),
    external_momenta(std::move(em_)),
    loc(lc_)
  {
  }


tree_expression& tree_expression::operator+=(const tree_expression& rhs)
  {
    if(!this->is_matching_type(rhs))
      throw exception(ERROR_COMPOSE_TREE_EXPRESSION_MISMATCHING_TYPE, exception_code::expression_error);

    // we know all variables and Wick product strings agree, so can just add the kernels
    this->K += rhs.K;

    return *this;
  }


void tree_expression::write(std::ostream& out) const
  {
    std::cout << "  time function = " << this->tm << '\n';
    std::cout << "  momentum kernel = " << this->K << '\n';
    std::cout << "  Wick product = " << this->WickProduct << '\n';
  }


bool tree_expression::is_matching_type(const tree_expression& obj) const
  {
    // test for equality of time function, Wick product, and external momenta
    const auto& at = this->tm;
    const auto& bt = obj.tm;

    if(!static_cast<bool>(at == bt)) return false;

    // test for equality of Wick product
    const auto& aw = this->WickProduct;
    const auto& bw = obj.WickProduct;

    if(!static_cast<bool>(aw == bw)) return false;

    // test for equality of external momenta
    auto a_em = order_symbol_set(this->external_momenta);
    auto b_em = order_symbol_set(obj.external_momenta);

    if(!std::equal(a_em.cbegin(), a_em.cend(), b_em.cbegin(), b_em.cend(),
                   [](const GiNaC::symbol& asym, const GiNaC::symbol& bsym) -> bool
                     { return asym.get_name() == bsym.get_name(); })) return false;

    return true;
  }


std::string tree_expression::to_Mathematica() const
  {
    std::ostringstream result;

    // write out time function and kernel
    result << "(" << format_print(this->tm) << ")*(" << this->K << ")";

    return result.str();
  }


tree_expression_key::tree_expression_key(const tree_expression& t)
  : tree(t)
  {
  }


size_t tree_expression_key::hash() const
  {
    // hash on time function, Wick product and external momenta

    // to hash the time expression, expand it completely and print
    // we are guaranteed that the expressions comes in a canonical order, even if that order is unpredictable
    std::ostringstream time_string;
    time_string << this->tree.get_time_function().expand();

    size_t h = 0;
    hash_impl::hash_combine(h, time_string.str());

    // to hash the Wick product, expand it completely and print
    std::ostringstream Wick_string;
    Wick_string << this->tree.get_Wick_product().expand();

    hash_impl::hash_combine(h, Wick_string.str());

    // order external momenta lexicographically, convert to a string, and hash
    auto ordered_em = order_symbol_set(this->tree.get_external_momenta());

    std::string em_string;
    std::for_each(ordered_em.begin(), ordered_em.end(),
                  [&](const GiNaC::symbol& e) -> void
                    { em_string += e.get_name(); });

    hash_impl::hash_combine(h, em_string);

    return h;
  }


bool tree_expression_key::is_equal(const tree_expression_key& obj) const
  {
    return this->tree.is_matching_type(obj.tree);
  }
