//
// Created by David Seery on 31/08/2017.
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


#include <map>

#include "legendre_utils.h"
#include "special_functions.h"

#include "utilities/GiNaC_utils.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


using CosToLegendreDatabase = std::map<unsigned int, GiNaC::matrix>;
using LegendrePolynomialDatabase = std::map<unsigned int, GiNaC::ex>;

//! database of cache mapping matrices
static CosToLegendreDatabase matrix_db;

//! database of Legendre polynomials
static LegendrePolynomialDatabase Pn_db;

//! symbol used internally to construct Legendre polynomials
static GiNaC::symbol x{"x"};


//! comppute Legendre polynomial of order n
const GiNaC::ex& LegendreP(unsigned int n)
  {
    auto t = Pn_db.find(n);
    if(t != Pn_db.end()) return t->second;

    if(n == 0) { Pn_db[0] = GiNaC::ex(1); return Pn_db[0]; }
    if(n == 1) { Pn_db[1] = x; return Pn_db[1]; }

    // use Bonnet recursion to generate higher polynomials
    GiNaC::numeric m(n);
    auto Pm = ((2*m-1)*x*LegendreP(n-1) - (m-1)*LegendreP(n-2)) / m;

    Pn_db[n] = Pm;
    return Pn_db[n];
  }


//! compute transformation matrix from powers of Cos to Legendre representation
const GiNaC::matrix& Cos_to_Legendre_matrix(unsigned int deg)
  {
    // do we have a cached copy of this matrix?
    auto t = matrix_db.find(deg);
    if(t != matrix_db.end()) return t->second;

    // no, so we have to compute it from scratch

    // populate a matrix with the coefficients of Pn as functions of x
    GiNaC::matrix coeffs{deg+1, deg+1};

    for(unsigned int i = 0; i <= deg; ++i)
      {
        // get ith Legendre polynomial
        const auto& Pn = LegendreP(i).expand();

        for(unsigned int j = 0; j <= deg; ++j)
          {
            coeffs(i, j) = Pn.coeff(x, j);
          }
      }

    // the required conversion matrix is the inverse of this coefficient matrix
    matrix_db[deg] = coeffs.inverse();

    return matrix_db[deg];
  }


//! compute Legendre polynomial of order n with argument t
GiNaC::ex LegP(unsigned int n, const GiNaC::ex& t)
  {
    const auto& Pn = LegendreP(n);
    GiNaC::exmap map;
    map[x] = t;

    return Pn.subs(map);
  }


namespace dot_products_to_cos_impl
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
            // does this factor carry an index? if not: convert it, push it to the result and carry on
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
          return GiNaC::pow(convert(base).expand(GiNaC::expand_options::expand_indexed), exponent);

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
        if(GiNaC::is_a<GiNaC::add>(expr))
          {
            return convert_add(expr);
          }
        if(GiNaC::is_a<GiNaC::mul>(expr))
          {
            return convert_mul(expr);
          }
        if(GiNaC::is_a<GiNaC::power>(expr))
          {
            return convert_power(expr);
          }

        // otherwise assume nothing to do, so return
        return expr;
      }

  }   // namespace dot_products_to_cos_impl


GiNaC::ex dot_products_to_cos(const GiNaC::ex& expr)
  {
    // step through the kernel, converting any explicit inner products into cosines;
    // these can later be converted into Legendre polynomials if required
    using dot_products_to_cos_impl::convert;

    return convert(expr.expand(GiNaC::expand_options::expand_indexed));
  }


namespace cosine_Legendre_impl
  {

    void find_pairs(const GiNaC::symbol& q, const GiNaC::ex& expr, GiNaC_symbol_set& set,
                    const std::string& name, unsigned int op1, unsigned int op2)
      {
        if(GiNaC::is_a<GiNaC::function>(expr))
          {
            const auto& f = GiNaC::ex_to<GiNaC::function>(expr);
            if(f.get_name() == name)
              {
                const auto& p1 = GiNaC::ex_to<GiNaC::symbol>(f.op(op1));
                const auto& p2 = GiNaC::ex_to<GiNaC::symbol>(f.op(op2));

                if(p1 == q) set.insert(p2);
                if(p2 == q) set.insert(p1);
                return;
              }
          }

        // arrive here if expr isn't a function, or the function wasn't what we were looking for
        for(size_t i = 0; i < expr.nops(); ++i)
          {
            find_pairs(q, expr.op(i), set, name, op1, op2);
          }
      }


    GiNaC_symbol_set get_Cos_pairs(const GiNaC::symbol& q, const GiNaC::ex& expr)
      {
        GiNaC_symbol_set pair_set;
        find_pairs(q, expr, pair_set, "Cos", 0, 1);

        return pair_set;
      }


    GiNaC_symbol_set get_LegP_pairs(const GiNaC::symbol& q, const GiNaC::ex& expr)
      {
        GiNaC_symbol_set pair_set;
        find_pairs(q, expr, pair_set, "LegP", 1, 2);

        return pair_set;
      }


    unsigned int get_max_LegP_order(const GiNaC::symbol& p1, const GiNaC::symbol& p2, const GiNaC::ex& expr)
      {
        if(GiNaC::is_a<GiNaC::function>(expr))
          {
            const auto& f = GiNaC::ex_to<GiNaC::function>(expr);
            if(f.get_name() == "LegP")
              {
                const auto& o = GiNaC::ex_to<GiNaC::numeric>(f.op(0));
                const auto& q1 = GiNaC::ex_to<GiNaC::symbol>(f.op(1));
                const auto& q2 = GiNaC::ex_to<GiNaC::symbol>(f.op(2));

                if(p1 == q1 && p2 == q2) return static_cast<unsigned int>(o.to_int());
                return 0;
              }
          }

        unsigned int res = 0;
        for(size_t i = 0; i < expr.nops(); ++i)
          {
            res = std::max(res, get_max_LegP_order(p1, p2, expr.op(i)));
          }

        return res;
      }


  }   // namespace cosine_Legendre_impl


GiNaC::ex cosines_to_Legendre(const GiNaC::ex& expr, const GiNaC::symbol& q)
  {
    using cosine_Legendre_impl::get_Cos_pairs;

    // obtain the set of symbols with which q appears in conjunction
    auto ex_exp = expr.expand();
    auto pair_set = get_Cos_pairs(q, ex_exp);

    // for each symbol, collect a polynomial in Cos(q,) and exchange it for a Legendre representation
    for(const auto& p : pair_set)
      {
        // arguments will be ordered canonically
        GiNaC::ex c;
        GiNaC::symbol p1 = q;
        GiNaC::symbol p2 = p;
        if(std::less<GiNaC::symbol>{}(p2, p1)) std::swap(p1, p2);
        c = Angular::Cos(p1, p2);

        auto deg = static_cast<unsigned int>(ex_exp.degree(c));
        if(deg == 0)
          throw exception(ERROR_COULDNT_COLLECT_COS, exception_code::loop_transformation_error);

        // collect coefficients of 1, cos, cos^2, ..., cos^deg into a row vector
        GiNaC::matrix coeffs{1, deg+1};

        for(unsigned int i = 0; i <= deg; ++i)
          {
            coeffs(0, i) = ex_exp.coeff(c, i);
          }

        // set up Cos-to-Legendre matrix of required degree
        const auto& xfm = Cos_to_Legendre_matrix(deg);

        // build column vector of Legendre polynomials
        GiNaC::matrix LegP{deg+1, 1};

        for(unsigned int i = 0; i <= deg; ++i)
          {
            LegP(i, 0) = Angular::LegP(GiNaC::numeric(i), p1, p2);
          }

        GiNaC::matrix res = coeffs.mul(xfm).mul(LegP);
        if(res.cols() != 1)
          throw exception(ERROR_LEGENDRE_TRANSFORM_UNEXPECTED_SIZE, exception_code::loop_transformation_error);
        if(res.cols() != 1)
          throw exception(ERROR_LEGENDRE_TRANSFORM_UNEXPECTED_SIZE, exception_code::loop_transformation_error);

        // store Legendre representation
        ex_exp = res(0,0).expand();
      }

    return ex_exp;
  }


GiNaC::ex Legendre_to_cosines(GiNaC::ex expr, const GiNaC::symbol q)
  {
    using cosine_Legendre_impl::get_LegP_pairs;
    using cosine_Legendre_impl::get_max_LegP_order;

    // obtain the set of symbols with which q appears in conjunction
    auto pair_set = get_LegP_pairs(q, expr);

    // for each symbol, substitute for the Legendre polynomials
    for(const auto& p : pair_set)
      {
        // arguments will be ordered canonically
        GiNaC::symbol p1 = q;
        GiNaC::symbol p2 = p;
        if(std::less<GiNaC::symbol>{}(p2, p1)) std::swap(p1, p2);

        // determine maximum order with which this combination occurs
        unsigned int max_order = get_max_LegP_order(p1, p2, expr);

        for(unsigned int i = 0; i <= max_order; i++)
          {
            GiNaC::exmap map;
            map[Angular::LegP(GiNaC::numeric(i), p1, p2)] = LegP(i, Angular::Cos(p1, p2));
            expr = expr.subs(map);
          }
      }

    return expr;
  }
