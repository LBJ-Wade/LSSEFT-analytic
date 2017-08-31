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
#include "detail/legendre_utils.h"
#include "detail/contractions.h"

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
    std::cout << "  Wick product = " << this->WickProduct << '\n';

    if(!this->loop_momenta.empty())
      {
        std::cout << "  loop momenta =";
        for(const auto& sym : this->loop_momenta)
          {
            std::cout << " " << sym;
          }
        std::cout << '\n';
      }

    if(!this->Rayleigh_momenta.empty())
      {
        std::cout << "  Rayleigh momenta =";
        for(const auto& rule : this->Rayleigh_momenta)
          {
            std::cout << " " << rule.first << " -> " << rule.second << ";";
          }
        std::cout << '\n';
      }
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


void loop_integral::dot_products_to_cos()
  {
    // step through the kernel, converting any explicit inner products into cosines;
    // these can later be converted into Legendre polynomials if required
    GiNaC::ex new_K = dot_products_to_cos_impl::convert(this->K.expand(GiNaC::expand_options::expand_indexed));
    this->K = new_K;
  }


void loop_integral::canonicalize_loop_labels()
  {
    GiNaC_symbol_set new_momenta;
    GiNaC::exmap relabel;

    // step through loop momenta, relabelling to canonicalized variables
    unsigned int count = 0;
    for(const auto& l : this->loop_momenta)
      {
        const auto L = this->sf.make_canonical_loop_momentum(count++);
        relabel[l] = L;
        new_momenta.insert(L);
      }

    // propagate this relabelling to the kernel and the Wick product
    this->K = this->K.subs(relabel);
    this->WickProduct = this->WickProduct.subs(relabel);

    // propagate to any Rayleigh replacement rules in use
    for(auto& rule : this->Rayleigh_momenta)
      {
        rule.second = rule.second.subs(relabel);
      }

    this->loop_momenta = new_momenta;
  }

void loop_integral::canonicalize_Rayleigh_labels()
  {
    GiNaC::exmap new_Rayleigh;
    GiNaC::exmap relabel;

    // step through Rayleigh momenta, relabelling to canonicalized forms
    unsigned int count = 0;
    for(const auto& rule : this->Rayleigh_momenta)
      {
        const auto S = this->sf.make_canonical_Rayleigh_momentum(count++);
        relabel[rule.first] = S;
        new_Rayleigh[S] = rule.second;    // RHS of rules never depend on other Rayleigh momenta, so this is safe
      }

    // propagate this relabelling to the kernel and Wick product
    this->K = this->K.subs(relabel);
    this->WickProduct = this->WickProduct.subs(relabel);

    this->Rayleigh_momenta = new_Rayleigh;
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


void loop_integral::cosines_to_Legendre(const GiNaC::symbol& q)
  {
    using cosine_Legendre_impl::get_Cos_pairs;
    
    // obtain the set of symbols with which q appears in conjunction
    auto pair_set = get_Cos_pairs(q, this->K);
    
    // for each symbol, collect a polynomial in Cos(q,) and exchange it for a Legendre representation
    for(const auto& p : pair_set)
      {
        // arguments will be ordered canonically
        GiNaC::ex c;
        GiNaC::symbol p1 = q;
        GiNaC::symbol p2 = p;
        if(std::less<GiNaC::symbol>{}(p2, p1)) std::swap(p1, p2);
        c = Angular::Cos(p1, p2);

        auto deg = static_cast<unsigned int>(this->K.degree(c));
        if(deg == 0)
          throw exception(ERROR_COULDNT_COLLECT_COS, exception_code::loop_transformation_error);
        
        // collect coefficients of 1, cos, cos^2, ..., cos^deg into a row vector
        GiNaC::matrix coeffs{1, deg+1};
        
        for(unsigned int i = 0; i <= deg; ++i)
          {
            coeffs(0, i) = this->K.coeff(c, i);
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
        this->K = res(0,0);
      }
  }


void loop_integral::Legendre_to_cosines(const GiNaC::symbol q)
  {
    using cosine_Legendre_impl::get_LegP_pairs;
    using cosine_Legendre_impl::get_max_LegP_order;

    // obtain the set of symbols with which q appears in conjunction
    auto pair_set = get_LegP_pairs(q, this->K);

    // for each symbol, substitute for the Legendre polynomials
    for(const auto& p : pair_set)
      {
        // arguments will be ordered canonically
        GiNaC::symbol p1 = q;
        GiNaC::symbol p2 = p;
        if(std::less<GiNaC::symbol>{}(p2, p1)) std::swap(p1, p2);

        // determine maximum order with which this combination occurs
        unsigned int max_order = get_max_LegP_order(p1, p2, this->K);

        for(unsigned int i = 0; i <= max_order; i++)
          {
            GiNaC::exmap map;
            map[Angular::LegP(GiNaC::numeric(i), p1, p2)] = LegP(i, Angular::Cos(p1, p2));
            this->K = this->K.subs(map);
          }
      }
  }


void loop_integral::reduce_angular_integrals()
  {
    // nothing to do at tree-level
    if(this->loop_momenta.empty()) return;

    // apply Rayleigh reduction algorithm at one-loop
    if(this->loop_momenta.size() == 1)
      {
        this->reduce_angular_integrals_one_loop();
        return;
      }
  }


void loop_integral::reduce_angular_integrals_one_loop()
  {
    // STEP 1 - canonicalize all labels
    this->canonicalize_loop_labels();
    this->canonicalize_Rayleigh_labels();

    // STEP 2 - convert all dot products to explicit cosines
    this->dot_products_to_cos();

    // STEP 3 - pair nontrivial arguments in Wick product with Rayleigh momenta
    this->match_Wick_to_Rayleigh();

    // STEP 4 - apply Rayleigh plane-wave expansion to the Rayleigh momentum
    if(this->Rayleigh_momenta.size() > 1)
      throw exception(ERROR_MULTIPLE_RAYLEIGH_MOMENTA_NOT_IMPLEMENTED, exception_code::loop_transformation_error);
  }


void loop_integral::match_Wick_to_Rayleigh()
  {
    GiNaC::ex new_Wick{1};

    if(!GiNaC::is_a<GiNaC::mul>(this->WickProduct))
      throw exception(ERROR_BADLY_FORMED_WICK_PRODUCT, exception_code::loop_transformation_error);

    for(size_t i = 0; i < this->WickProduct.nops(); ++i)
      {
        const auto& factor = this->WickProduct.op(i);
        const auto& f = GiNaC::ex_to<GiNaC::function>(factor);

        if(f.get_name() == "Pk")
          {
            const auto& f1 = GiNaC::ex_to<GiNaC::symbol>(f.op(0));
            const auto& f2 = GiNaC::ex_to<GiNaC::symbol>(f.op(1));
            const auto& arg = f.op(2);

            // if argument is a simple symbol (a loop momentum or a Rayleigh momentum, or an external momentum) then we have nothing to do
            if(GiNaC::is_a<GiNaC::symbol>(arg))
              {
                const auto& sym = GiNaC::ex_to<GiNaC::symbol>(arg);
                if(this->loop_momenta.find(sym) != this->loop_momenta.end()) { new_Wick *= cfs::Pk(f1, f2, sym); continue; }
                if(this->Rayleigh_momenta.find(sym) != this->Rayleigh_momenta.end()) { new_Wick *= cfs::Pk(f1, f2, sym); continue; }
                if(this->external_momenta.find(sym) != this->external_momenta.end()) { new_Wick *= cfs::Pk(f1, f2, sym); continue; }

                std::cerr << factor << '\n';

                throw exception(ERROR_UNKNOWN_WICK_PRODUCT_LABEL, exception_code::loop_transformation_error);
              }

            // otherwise, try to match this argument to a Rayleigh momentum (or its negative)
            auto t = this->Rayleigh_momenta.begin();
            for(; t != this->Rayleigh_momenta.end(); ++t)
              {
                // if a direct match, just relabel the argument; don't need to worry about sign inversion
                // since the power spectrum just depends on |arg|
                if(static_cast<bool>(t->second == arg)
                   || static_cast<bool>(t->second == -arg))
                  {
                    new_Wick *= cfs::Pk(f1, f2, t->first);
                    break;
                  }
              }

            if(t == this->Rayleigh_momenta.end())
              {
                if(!this->Rayleigh_momenta.empty())
                  {
                    std::cerr << factor << '\n';
                    throw exception(ERROR_CANT_MATCH_WICK_TO_RAYLEIGH, exception_code::loop_transformation_error);
                  }

                // matching Rayleigh momentum didn't exist, so insert one
                auto S = this->sf.make_canonical_Rayleigh_momentum(0);
                new_Wick *= cfs::Pk(f1, f2, S);
                this->Rayleigh_momenta[S] = arg;
              }

            continue;
          }

        throw exception(ERROR_BADLY_FORMED_WICK_PRODUCT, exception_code::loop_transformation_error);
      }

    // replace old Wick product with matched version
    this->WickProduct = new_Wick;
  }
