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

#ifndef LSSEFT_ANALYTIC_TREE_DB_H
#define LSSEFT_ANALYTIC_TREE_DB_H


#include <iostream>
#include <unordered_map>
#include <memory>

#include "lib/expressions/tree_expression.h"


//! a tree_db is a container for the 0-loop integrals generated as part of a correlation function computation
class tree_db
  {

    // TYPES

  public:

    //! expr_type
    using expr_type = tree_expression;

  protected:

    //! database is a set of tree_expressions, keyed by tree_expression_key
    using db_type = std::unordered_map< tree_expression_key, std::unique_ptr<tree_expression> >;

  public:

    //! iterator
    using iterator = db_type::iterator;

    //! const iterator
    using const_iterator = db_type::const_iterator;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor is default
    tree_db() = default;

    //! copy constructor
    tree_db(const tree_db& obj);

    //! destructor is default
    ~tree_db() = default;


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
    void emplace(std::unique_ptr<tree_expression> elt);


    // ARITHMETIC

  public:

    //! increment
    tree_db& operator+=(const tree_db& obj);


    // TRANSFORMATIONS

  public:

    //! apply simplification map
    void simplify(const GiNaC::exmap& map);

    //! canonicalize external momenta by converting any angular factors involving them to cosines
    //! rather than Legendre polynomials
    void canonicalize_external_momenta();


    // SERVICES

  public:

    //! write self to steam
    void write(std::ostream& out) const;

    //! write Mathematica-format expression
    void write_Mathematica(std::ostream& out, std::string symbol) const;


    // INTERNAL DATA

  protected:

    //! database of tree-level expressions
    db_type db;

  };


//! perform stream insertion
std::ostream& operator<<(std::ostream& str, const tree_db& obj);


#endif //LSSEFT_ANALYTIC_TREE_DB_H
