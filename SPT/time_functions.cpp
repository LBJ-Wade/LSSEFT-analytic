//
// Created by David Seery on 23/08/2017.
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

#include "time_functions.h"


GiNaC::ex D_deriv(const GiNaC::ex& z, unsigned int diff_param)
  {
    return -SPT::D(z)*SPT::f(z) / (1+z);
  }


GiNaC::ex DA_deriv(const GiNaC::ex& z, unsigned int diff_param)
  {
    return -SPT::DA(z)*SPT::fA(z) / (1+z);
  }


GiNaC::ex DB_deriv(const GiNaC::ex& z, unsigned int diff_param)
  {
    return -SPT::DB(z)*SPT::fB(z) / (1+z);
  }


GiNaC::ex DD_deriv(const GiNaC::ex& z, unsigned int diff_param)
  {
    return -SPT::DD(z)*SPT::fD(z) / (1+z);
  }


GiNaC::ex DE_deriv(const GiNaC::ex& z, unsigned int diff_param)
  {
    return -SPT::DE(z)*SPT::fE(z) / (1+z);
  }


GiNaC::ex DF_deriv(const GiNaC::ex& z, unsigned int diff_param)
  {
    return -SPT::DF(z)*SPT::fF(z) / (1+z);
  }


GiNaC::ex DG_deriv(const GiNaC::ex& z, unsigned int diff_param)
  {
    return -SPT::DG(z)*SPT::fG(z) / (1+z);
  }


GiNaC::ex DJ_deriv(const GiNaC::ex& z, unsigned int diff_param)
  {
    return -SPT::DJ(z)*SPT::fJ(z) / (1+z);
  }


namespace FRW
  {

    REGISTER_FUNCTION(Hub, dummy())
    
  }   // namespace FRW


namespace SPT
  {

    REGISTER_FUNCTION(f, dummy())
    REGISTER_FUNCTION(fA, dummy())
    REGISTER_FUNCTION(fB, dummy())
    REGISTER_FUNCTION(fD, dummy())
    REGISTER_FUNCTION(fE, dummy())
    REGISTER_FUNCTION(fF, dummy())
    REGISTER_FUNCTION(fG, dummy())
    REGISTER_FUNCTION(fJ, dummy())
    
    REGISTER_FUNCTION(D, derivative_func(D_deriv));
    REGISTER_FUNCTION(DA, derivative_func(DA_deriv));
    REGISTER_FUNCTION(DB, derivative_func(DB_deriv));
    REGISTER_FUNCTION(DD, derivative_func(DD_deriv));
    REGISTER_FUNCTION(DE, derivative_func(DE_deriv));
    REGISTER_FUNCTION(DF, derivative_func(DF_deriv));
    REGISTER_FUNCTION(DG, derivative_func(DG_deriv));
    REGISTER_FUNCTION(DJ, derivative_func(DJ_deriv));
    
  }   // namespace SPT