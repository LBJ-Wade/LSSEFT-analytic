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

#include "special_functions.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"

namespace special
  {

    REGISTER_FUNCTION(j, dummy());

  }   // namespace special


namespace Angular
  {

    REGISTER_FUNCTION(Cos, dummy());
    REGISTER_FUNCTION(LegP, dummy());

  }   // namespace Angular


namespace Fabrikant
  {

    GiNaC::ex FabJ_eval(const GiNaC::ex& lambda, const GiNaC::ex& mu, const GiNaC::ex& nu,
                        const GiNaC::ex& s, const GiNaC::ex& t, const GiNaC::ex& u)
      {
        const auto& lambda_num = GiNaC::ex_to<GiNaC::numeric>(lambda);
        const auto& mu_num = GiNaC::ex_to<GiNaC::numeric>(mu);
        const auto& nu_num = GiNaC::ex_to<GiNaC::numeric>(nu);

        if(lambda_num.to_int() != 0)
          throw exception(ERROR_FABJ_FIRST_ARG_NONZERO, exception_code::Fabrikant_error);

        if(mu_num.to_int() != nu_num.to_int())
          throw exception(ERROR_FABJ_SECOND_ARGS_NOT_EQUAL, exception_code::Fabrikant_error);

        if(mu_num.to_int() < 0)
          throw exception(ERROR_FABJ_SECOND_ARGS_NEGATIVE, exception_code::Fabrikant_error);

        if(mu_num.to_int() == 0)
          return GiNaC::Pi/(4*s*t*u);

        if(mu_num.to_int() == 1)
          return (GiNaC::Pi*(-GiNaC::pow(s,2) + GiNaC::pow(t,2) + GiNaC::pow(u,2)))/(8*s*GiNaC::pow(t,2)*GiNaC::pow(u,2));

        if(mu_num.to_int() == 2)
          return (GiNaC::Pi*(3*GiNaC::pow(GiNaC::pow(s,2) - GiNaC::pow(t,2),2) + 2*(-3*GiNaC::pow(s,2) + GiNaC::pow(t,2))*GiNaC::pow(u,2) +
                      3*GiNaC::pow(u,4)))/(32*s*GiNaC::pow(t,3)*GiNaC::pow(u,3));

        return FabJ(lambda, mu, nu, s, t, u).hold();
      }

    REGISTER_FUNCTION(FabJ, eval_func(FabJ_eval));

  }   // namespace Fabrikant
