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

#include "loop_integral.h"

#include "shared/common.h"
#include "shared/exceptions.h"
#include "utilities/GiNaC_utils.h"
#include "utilities/hash_combine.h"
#include "localizations/messages.h"


namespace one_loop_reduced_integral_impl
  {

    //! forward-declare class key
    class key;

    //! integration element captures a single integration
    class integration_element
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor captures integrand, measure, integration variables, Wick product, time factor, external momenta
        integration_element(GiNaC::ex ig_, GiNaC::ex ms_, GiNaC::ex wp_, time_function tm_, GiNaC_symbol_set vs_,
                            GiNaC_symbol_set em_);

        //! destructor is default
        ~integration_element() = default;


        // ACCESSORS

      public:

        //! get integration variable list
        const GiNaC_symbol_set& get_integration_variables() const { return this->variables; }


        // SERVICES

      public:

        //! write self to stream
        void write(std::ostream& str) const;


        // TRANSFORMATIONS

      public:

        //! apply simplification map
        void simplify(const GiNaC::exmap& map);

        //! convert external momenta to canonical Cos form
        void canonicalize_external_momenta();


        // UV LIMIT

      public:

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

        //! set of integration variabkes
        GiNaC_symbol_set variables;

        //! set of external momenta
        GiNaC_symbol_set external_momenta;


        friend class key;

      };


    //! perform stream insertion
    std::ostream& operator<<(std::ostream& str, const integration_element& obj);


    //! key is a flyweight class that indexes integration elements by their
    //! integration variables
    class key
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor accepts a reference to a integration_element
        explicit key(const integration_element& elt_);

        //! destructor is default
        ~key() = default;


        // SERVICES

      public:

        //! hash
        size_t hash() const;

        //! compare for equality
        bool is_equal(const key& obj) const;


        // INTERNAL DATA

      private:

        //! cache reference to partner class
        const integration_element& elt;

      };


    //! an integrand database is a map from integration variables to a list of integrands
    using integrand_db = std::unordered_map< key, std::vector< std::unique_ptr<integration_element> > >;

  }   // namespace one_loop_reduced_integral


// specialize std::hash<> and std::is_equal<> to key
namespace std
  {

    template <>
    struct hash<one_loop_reduced_integral_impl::key>
      {
        size_t operator()(const one_loop_reduced_integral_impl::key& obj) const
          {
            return obj.hash();
          }
      };

    template <>
    struct equal_to<one_loop_reduced_integral_impl::key>
      {
        bool operator()(const one_loop_reduced_integral_impl::key& a, const one_loop_reduced_integral_impl::key& b) const
          {
            return a.is_equal(b);
          }
      };

  }   // namespace std


class one_loop_reduced_integral
  {

    // TYPES

  protected:

    //! pull in integrand_db
    using integrand_db = one_loop_reduced_integral_impl::integrand_db;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor accepts a loop_integral and performs dimensional reduction on it
    one_loop_reduced_integral(const loop_integral& i_, symbol_factory& sf_);

    //! destructor is default
    ~one_loop_reduced_integral() = default;


    // TRANSFORMATIONS

  public:

    //! apply simplification map
    void simplify(const GiNaC::exmap& map);

    //! canonicalize external momenta
    void canonicalize_external_momenta();

  protected:

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


    // UV LIMIT

  public:

    //! get UV limit
    GiNaC::ex get_UV_limit(unsigned int order=2) const;


    // SERVICES

  public:

    //! write self to stream
    void write(std::ostream& out) const;


    // INTERNAL DATA

  private:

    // AGENTS

    //! cache reference to symbol factory
    symbol_factory& sf;


    // DATA

    //! cache reference to original loop integral
    const loop_integral& loop_int;

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


    // SPECIAL SYMBOLS

    //! integration variable x
    const GiNaC::symbol x;


    // REDUCED EXPRESSIONS

    //! set up an integrand database
    integrand_db integrand;

  };


template <typename VisitorFunction>
GiNaC::ex one_loop_reduced_integral::integrate_Legendre(const GiNaC::ex& term, const GiNaC::symbol& q, VisitorFunction f)
  {
    auto term_ex = term.expand();

    if(GiNaC::is_a<GiNaC::mul>(term_ex))
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

    throw exception(ERROR_BADLY_FORMED_TOP_LEVEL_LEGENDRE_SUM, exception_code::loop_transformation_error);
  }


//! perform stream insertion
std::ostream& operator<<(std::ostream& str, const one_loop_reduced_integral& obj);


#endif //LSSEFT_ANALYTIC_ONE_LOOP_REDUCED_INTEGRAL_H
