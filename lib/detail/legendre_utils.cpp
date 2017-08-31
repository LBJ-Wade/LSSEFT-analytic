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

#include "legendre_utils.h"

#include <map>

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
        const auto& Pn = LegendreP(i);

        for(unsigned int j = 0; j <= deg; ++j)
          {
            coeffs(i, j) = Pn.coeff(x, j);
          }
      }

    // the required conversion matrix is the inverse of this coefficient matrix
    matrix_db[deg] = coeffs.inverse();

    return matrix_db[deg];
  }

GiNaC::ex LegP(unsigned int n, const GiNaC::ex& t)
  {
    const auto& Pn = LegendreP(n);
    GiNaC::exmap map;
    map[x] = t;

    return Pn.subs(map);
  }
