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

#ifndef LSSEFT_ANALYTIC_ONE_LOOP_REDUCED_INTEGRAL_H
#define LSSEFT_ANALYTIC_ONE_LOOP_REDUCED_INTEGRAL_H


#include <unordered_map>

#include "lib/expressions/oneloop_expression.h"

#include "shared/common.h"
#include "shared/exceptions.h"
#include "utilities/GiNaC_utils.h"
#include "utilities/hash_combine.h"
#include "localizations/messages.h"


//! forward-declare class key
class oneloop_element_key;



//! integration element captures a single integration
class oneloop_element
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor captures integrand, measure, integration variables, Wick product, time factor, external momenta
    oneloop_element(GiNaC::ex ig_, GiNaC::ex ms_, GiNaC::ex wp_, time_function tm_,
                    GiNaC_symbol_set vs_, GiNaC::symbol ang_, GiNaC_symbol_set em_);

    //! destructor is default
    ~oneloop_element() = default;


    // OPERATIONS

  public:

    //! increment in-place
    oneloop_element& operator+=(const oneloop_element& rhs);


    // ACCESSORS

  public:

    //! get integrand
    const GiNaC::ex& get_integrand() const { return this->integrand; }

    //! get measure
    const GiNaC::ex& get_measure() const { return this->measure; }

    //! get Wick product
    const GiNaC::ex& get_Wick_product() const { return this->WickProduct; }

    //! get time function
    const GiNaC::ex& get_time_function() const { return this->tm; }

    //! get integration variable list
    const GiNaC_symbol_set& get_integration_variables() const { return this->variables; }

    //! get external momenta
    const GiNaC_symbol_set& get_external_momenta() const { return this->external_momenta; }


    // SERVICES

  public:

    //! write self to stream
    void write(std::ostream& str) const;

    //! test for matching type (matches time function, measure, Wick product, integration variables, external momenta)
    bool is_matching_type(const oneloop_element& obj) const;

    //! test for nullity of integrand
    bool null() const { return static_cast<bool>(this->integrand == 0); }


    // TRANSFORMATIONS

  public:

    //! apply simplification map
    void simplify(const GiNaC::exmap& map);

    //! convert external momenta to canonical Cos form
    void canonicalize_external_momenta();

    //! filter integrand
    void filter(const GiNaC::symbol& pattern, unsigned int order = 1);


    // SERVICES

  public:

    //! convert self to Mathematica format expression
    std::string to_Mathematica(bool do_dx) const;

    //! construct UV limit
    GiNaC::ex get_UV_limit(unsigned int order=2) const;


    // INTERNAL DATA

  private:

    //! integrand
    GiNaC::ex integrand;

    //! measure
    GiNaC::ex measure;

    //! Wick product
    GiNaC::ex WickProduct;

    //! time function
    time_function tm;

    //! set of integration variables
    GiNaC_symbol_set variables;

    //! set of external momenta
    GiNaC_symbol_set external_momenta;


    //! cache angular integration variable
    const GiNaC::symbol angular_dx;


    friend class oneloop_element_key;

  };


//! perform stream insertion
std::ostream& operator<<(std::ostream& str, const oneloop_element& obj);


//! key is a flyweight class that indexes integration elements by their
//! integration variables
class oneloop_element_key
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor accepts a reference to a integration_element
    explicit oneloop_element_key(const oneloop_element& elt_);

    //! destructor is default
    ~oneloop_element_key() = default;


    // SERVICES

  public:

    //! return time function associated with this key
    const time_function& get_time_function() const;

    //! hash
    size_t hash() const;

    //! compare for equality
    bool is_equal(const oneloop_element_key& obj) const;


    // INTERNAL DATA

  private:

    //! cache reference to partner class
    const oneloop_element& elt;

  };


//! an integrand database is a map from integration with given dummy variables, time function and Wick product
//! to the cumulative integrand
using oneloop_element_db = std::unordered_map< oneloop_element_key, std::unique_ptr<oneloop_element> >;

//! stream insertion
std::ostream& operator<<(std::ostream& out, const oneloop_element_db& obj);


// specialize std::hash<> and std::is_equal<> to key
namespace std
  {

    template <>
    struct hash<oneloop_element_key>
      {
        size_t operator()(const oneloop_element_key& obj) const
          {
            return obj.hash();
          }
      };

    template <>
    struct equal_to<oneloop_element_key>
      {
        bool operator()(const oneloop_element_key& a, const oneloop_element_key& b) const
          {
            return a.is_equal(b);
          }
      };

  }   // namespace std


class oneloop_reduced_integral
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor accepts a loop_integral container and performs dimensional reduction on it
    oneloop_reduced_integral(const oneloop_expression& i_, service_locator& lc_, bool s_);

    //! destructor is default
    ~oneloop_reduced_integral() = default;


    // ACCESSORS

  public:

    //! get database
    const oneloop_element_db& get_db() const { return this->integrand; }

    //! query empty status
    bool empty() const { return this->integrand.empty(); }


    // TRANSFORMATIONS

  public:

    //! apply simplification map
    void simplify(const GiNaC::exmap& map);

    //! canonicalize external momenta
    void canonicalize_external_momenta();

    //! remove empty records
    void prune();

  protected:

    //! emplace a one_loop_element in the database
    void emplace(std::unique_ptr<oneloop_element> elt);

    //! apply one-loop reduction formula to a single term
    void reduce(const GiNaC::ex& expr);

    //! apply one-loop reduction to a term with zero Rayleigh momenta
    void one_loop_reduce_zero_Rayleigh(const GiNaC::ex& term);

    //! apply one-loop reduction to a term with one Rayleigh momentum
    void one_loop_reduce_one_Rayleigh(const GiNaC::ex& term, const GiNaC::symbol& R);

    //! generic function to walk an expanded expression, applying a function to integrate Legendre products
    template <typename VisitorFunction>
    GiNaC::ex integrate_Legendre(const GiNaC::ex& term, const GiNaC::symbol& q, VisitorFunction f);

    //! apply Legendre orthogonality formula to integrate a term containing zero, one or two Legendre polynomials
    GiNaC::ex apply_Legendre_orthogonality(const GiNaC::ex& expr, const GiNaC::symbol& q);

    //! apply Legendre orthogonality formula to a one-loop term generated by Rayleigh momentum.
    //! L is the loop momentum and k is an external momentum
    GiNaC::ex apply_Legendre_orthogonality(const GiNaC::ex& expr, const GiNaC::symbol& L, const GiNaC::numeric& Lcoeff,
                                           const GiNaC::symbol& k, const GiNaC::numeric& kcoeff, const GiNaC::symbol& R);


    // SERVICES

  public:

    //! write self to stream
    void write(std::ostream& out) const;

    //! format self as Mathematica expressions
    std::string to_Mathematica(bool do_dx) const;

    //! get UV limit
    GiNaC::ex get_UV_limit(unsigned int order=2) const;


    // INTERNAL DATA

  private:

    // AGENTS

    //! cache reference to service locator
    service_locator& loc;


    // DATA

    //! cache reference to original loop integral
    const oneloop_expression& loop_int;

    //! cache list of Rayleigh momenta
    const subs_list& Rayleigh_momenta;

    //! cache Wick product
    const GiNaC::ex& WickProduct;

    //! cache time function
    const time_function& tm;

    //! cache loop momentum
    GiNaC::symbol loop_q;

    //! cache external momenta
    const GiNaC_symbol_set& external_momenta;

    //! symmetrize expession? - used for 22-type integrals
    bool symmetrize;


    // SPECIAL SYMBOLS

    //! integration variable x
    const GiNaC::symbol x;


    // REDUCED EXPRESSIONS

    //! set up an integrand database
    oneloop_element_db integrand;

  };


template <typename VisitorFunction>
GiNaC::ex oneloop_reduced_integral::integrate_Legendre(const GiNaC::ex& term, const GiNaC::symbol& q, VisitorFunction f)
  {
    auto term_ex = term.expand();

    if(GiNaC::is_a<GiNaC::mul>(term_ex) || GiNaC::is_a<GiNaC::numeric>(term_ex) || GiNaC::is_a<GiNaC::power>(term_ex))
      {
        return f(term_ex, q);
      }

    if(GiNaC::is_a<GiNaC::add>(term_ex))
      {
        GiNaC::ex temp{0};
        for(size_t i = 0; i < term_ex.nops(); ++i)
          {
            temp += f(term_ex.op(i), q);
          }
        return temp;
      }

    if(GiNaC::is_a<GiNaC::symbol>(term_ex) || GiNaC::is_a<GiNaC::function>(term_ex)) return term;

    std::cerr << term_ex << '\n';
    throw exception(ERROR_BADLY_FORMED_TOP_LEVEL_LEGENDRE_SUM, exception_code::loop_transformation_error);
  }


//! perform stream insertion
std::ostream& operator<<(std::ostream& out, const oneloop_reduced_integral& obj);


#endif //LSSEFT_ANALYTIC_ONE_LOOP_REDUCED_INTEGRAL_H
