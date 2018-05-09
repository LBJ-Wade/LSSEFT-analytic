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

#ifndef LSSEFT_ANALYTIC_TIME_FUNCTIONS_H
#define LSSEFT_ANALYTIC_TIME_FUNCTIONS_H


#include "ginac/ginac.h"


//! compute dD/dz
GiNaC::ex D_deriv(const GiNaC::ex& expr, unsigned int diff_param);

//! compute dDA/dz
GiNaC::ex DA_deriv(const GiNaC::ex& expr, unsigned int diff_param);

//! compute dDB/dz
GiNaC::ex DB_deriv(const GiNaC::ex& expr, unsigned int diff_param);

//! compute dDD/dz
GiNaC::ex DD_deriv(const GiNaC::ex& expr, unsigned int diff_param);

//! compute dDE/dz
GiNaC::ex DE_deriv(const GiNaC::ex& expr, unsigned int diff_param);

//! compute dDF/dz
GiNaC::ex DF_deriv(const GiNaC::ex& expr, unsigned int diff_param);

//! compute dDG/dz
GiNaC::ex DG_deriv(const GiNaC::ex& expr, unsigned int diff_param);

//! compute dDJ/dz
GiNaC::ex DJ_deriv(const GiNaC::ex& expr, unsigned int diff_param);


// declare background cosmology functions
namespace FRW
  {
    
    //! Hubble rate H
    DECLARE_FUNCTION_1P(Hub)
    
  }

// declare time-dependent functions needed for 1-loop SPT

namespace SPT
  {
    
    //! linear growth factor f
    DECLARE_FUNCTION_1P(f)
    
    //! 1-loop growth function fA
    DECLARE_FUNCTION_1P(fA)
    
    //! 1-loop growth function fB
    DECLARE_FUNCTION_1P(fB)
    
    //! 1-loop growth function fD
    DECLARE_FUNCTION_1P(fD)
    
    //! 1-loop growth function fE
    DECLARE_FUNCTION_1P(fE)
    
    //! 1-loop growth function fF
    DECLARE_FUNCTION_1P(fF)
    
    //! 1-loop growth function fG
    DECLARE_FUNCTION_1P(fG)
    
    //! 1-loop growth function fJ
    DECLARE_FUNCTION_1P(fJ)

    
    //! linear growth function D
    DECLARE_FUNCTION_1P(D)
    
    //! 1-loop growth function DA
    DECLARE_FUNCTION_1P(DA);
    
    //! 1-loop growth function DB
    DECLARE_FUNCTION_1P(DB);
    
    //! 1-loop growth function DD
    DECLARE_FUNCTION_1P(DD);
    
    //! 1-loop growth function DE
    DECLARE_FUNCTION_1P(DE);
    
    //! 1-loop growth function DF
    DECLARE_FUNCTION_1P(DF);
    
    //! 1-loop growth function DG
    DECLARE_FUNCTION_1P(DG);
    
    //! 1-loop growth function DJ
    DECLARE_FUNCTION_1P(DJ);
    
  }   // namespace SPT



#endif //LSSEFT_ANALYTIC_TIME_FUNCTIONS_H
