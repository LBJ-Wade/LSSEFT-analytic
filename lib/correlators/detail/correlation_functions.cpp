//
// Created by David Seery on 06/04/2018.
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

#include "correlation_functions.h"


namespace cfs
  {

    //! Taylor expansion of Pk needs special attention, because we don't want to end up doing
    //! an expansion of something Pk(k).
    //! On the other hand, Pk(sqrt(k^2 + K^2 - 2kLx)) can certainly be expanded
    GiNaC::ex Pk_series(const GiNaC::ex& a, const GiNaC::ex& b, const GiNaC::ex& x,
                        const GiNaC::relational& rel, int order, unsigned options);

    // register Pk with custom series expansion function
    REGISTER_FUNCTION(Pk, series_func(Pk_series))

    GiNaC::ex Pk_series(const GiNaC::ex& a, const GiNaC::ex& b, const GiNaC::ex& x,
                        const GiNaC::relational& rel, int order, unsigned options)
      {
        // get symbol that we're expanding with respect to
        const auto& sym = GiNaC::ex_to<GiNaC::symbol>(rel.lhs());

        // if series expansion of argument has a constant term, then it's safe to do a normal Taylor
        // expansion
        if(x.coeff(sym, 0) != 0) throw GiNaC::do_taylor();

        // otherwise, we should avoid it because we will end up with an expansion of Pk(x) around x=0
        // instead, just return the original Pk(x) as the zeroth term of a power series
        return GiNaC::pseries(rel, GiNaC::epvector{ {Pk(a, b, x), 0} });
      }

  }   // namespace cfs = correlation functions
