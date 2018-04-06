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


#ifndef LSSEFT_ANALYTIC_ONELOOP_DB_H
#define LSSEFT_ANALYTIC_ONELOOP_DB_H


#include <iostream>
#include <unordered_map>
#include <memory>

#include "lib/expressions/oneloop_expression.h"
#include "lib/expressions/oneloop_reduced_integral.h"


//! a oneloop_db is a container for the 1-loop integrals generated as part of a correlation function computation
class oneloop_db
  {

    // TYPES

  public:

    //! a oneloop_pair is a doublet that groups a raw integral with its reduced form
    using oneloop_pair = std::pair< std::unique_ptr<oneloop_expression>, std::unique_ptr<oneloop_reduced_integral> >;

    //! expr_type
    using expr_type = oneloop_expression;

    //! reduced_type
    using reduced_type = oneloop_reduced_integral;

    //! pointer to reduced type
    using reduced_ptr_type = std::unique_ptr<reduced_type>;

  protected:

    //! database is a set of loop_pairs, keyed by oneloop_expression_key
    using db_type = std::unordered_map< oneloop_expression_key, oneloop_pair >;

  public:

    //! iterator
    using iterator = db_type::iterator;

    //! const iterator
    using const_iterator = db_type::const_iterator;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor is default
    oneloop_db() = default;

    //! copy constructor
    oneloop_db(const oneloop_db& obj);

    //! destructor is default
    ~oneloop_db() = default;


    // ITERATORS

  public:

    //! iterators
    iterator       begin()        { return this->db.begin(); }
    iterator       end()          { return this->db.end(); }
    const_iterator begin()  const { return this->db.cbegin(); }
    const_iterator end()    const { return this->db.cend(); }

    const_iterator cbegin() const { return this->db.cbegin(); }
    const_iterator cend()   const { return this->db.cend(); }


    // EMPLACE AN ELEMENT

  public:

    //! emplace an element
    void emplace(std::unique_ptr<oneloop_expression> elt);


    // FLUSH REDUCED INTEGRALS

  public:

    //! clear all reduced integrals
    void clear_reduced_integrals();


    // ARITHMETIC

  public:

    //! increment
    oneloop_db& operator+=(const oneloop_db& obj);


    // TRANSFORMATIONS

  public:

    //! reduce angular integrals
    void reduce_angular_integrals(service_locator& loc, bool symmetrize);

    //! apply simplification map
    void simplify(const GiNaC::exmap& map);

    //! canonicalize external momenta by converting any angular factors involving them to cosines
    //! rather than Legendre polynomials
    void canonicalize_external_momenta();

    //! prune empty records
    void prune();


    // SERVICES

  public:

    //! write self to steam
    void write(std::ostream& out) const;

    //! write Mathematica-format expression
    void write_Mathematica(std::ostream& out, std::string symbol, bool do_dx) const;

    //! compute UV limit
    GiNaC::ex get_UV_limit(unsigned int order=2) const;


    // INTERNAL DATA

  protected:

    //! database of 1-loop expressions
    db_type db;

  };


//! perform stream insertion
std::ostream& operator<<(std::ostream& str, const oneloop_db& obj);


#endif //LSSEFT_ANALYTIC_ONELOOP_DB_H
