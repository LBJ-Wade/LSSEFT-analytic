//
// Created by David Seery on 02/10/2017.
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


#include <sstream>

#include "LSSEFT.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


LSSEFT& LSSEFT::add(const Pk_rsd& P, std::string name)
  {
    auto it = this->db.find(name);

    if(it != this->db.end())
      {
        std::ostringstream msg;
        msg << ERROR_BACKEND_PK_RSD_ALREADY_REGISTERED_A
            << " '" << name << "' "
            << ERROR_BACKEND_PK_RSD_ALREADY_REGISTERED_B;
        throw exception(msg.str(), exception_code::backend_error);
      }

    auto res = this->db.insert(std::make_pair(std::move(name), std::cref(P)));
    if(!res.second) throw exception(ERROR_BACKEND_PK_INSERT_FAILED, exception_code::backend_error);

    return *this;
  }
