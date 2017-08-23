//
// Created by David Seery on 22/08/2017.
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
#ifndef LSSEFT_ANALYTIC_VECTOR_H
#define LSSEFT_ANALYTIC_VECTOR_H


#include "services/symbol_factory.h"

#include "ginac/ginac.h"


// forward-declare vector class
class vector;


// forward declare out-of-class operator overloads

//! unary - on a vector
vector operator-(const vector& a);

//! add two vectors
vector operator+(const vector& a, const vector& b);

//! subtract two vectors
vector operator-(const vector& a, const vector& b);

//! multiply a vector by a numeric: vector = LHS
vector operator*(const vector& a, const GiNaC::numeric b);

//! multiply a vector by a numeric: vector = RHS
vector operator*(const GiNaC::numeric a, const vector& b);

//! divide a vector by a numeric
vector operator/(const vector& a, const GiNaC::numeric b);

//! generate dot-product of two vectors
GiNaC::ex dot(const vector a, const vector b);


class vector
  {
  
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! no default public constructor; must be constructed through symbol_factory
    vector() = delete;
    
    //! destructor is default
    ~vector() = default;
    
  protected:
    
    //! basic private constructor accepts a symbol representing the vector kernel
    vector(GiNaC::symbol& k, symbol_factory& sf_);
    
    //! private constructor available to overloaded arithmetic functions that can accept
    //! a GiNaC::ex rather than just a symbol
    vector(GiNaC::ex e, symbol_factory& sf_);
    
    
    // INTERFACE
    
  public:
    
    //! get raw expression
    const GiNaC::ex& get_expr() const { return this->expr; }
    
    //! get component of a vector
    GiNaC::ex operator[](const GiNaC::idx& i) const;
    
    //! get norm-square of vector
    GiNaC::ex norm_square() const;
    
    //! get norm of vector
    GiNaC::ex norm() const;
    
    
    // INTERNAL DATA
    
  private:
    
    //! GiNaC expression representing a vector
    //! (not that this can be a linear combination of fundamental vectors)
    GiNaC::ex expr;
    
    //! cache reference to symbol factory
    symbol_factory& sf;
    
    
    friend class symbol_factory;
    
    friend vector operator-(const vector& a);
    friend vector operator+(const vector& a, const vector& b);
    friend vector operator-(const vector& a, const vector& b);
    friend vector operator*(const vector& a, const GiNaC::numeric b);
    friend vector operator*(const GiNaC::numeric a, const vector& b);
    friend vector operator/(const vector& a, const GiNaC::numeric b);
    friend GiNaC::ex dot(const vector a, const vector b);
    
  };


#endif //LSSEFT_ANALYTIC_VECTOR_H
