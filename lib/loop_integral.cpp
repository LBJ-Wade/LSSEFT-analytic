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

#include "loop_integral.h"


std::ostream& operator<<(std::ostream& str, const loop_integral& obj)
  {
    obj.write(str);
    return str;
  }


void loop_integral::write(std::ostream& out) const
  {
    std::cout << "  time function = " << this->tm << '\n';
    std::cout << "  momentum kernel = " << this->K << '\n';
    std::cout << "  Wick product = " << this->WickString << '\n';

    std::cout << "  loop momenta =";
    for(const auto& sym : this->loop_momenta)
      {
        std::cout << " " << sym;
      }
    if(this->loop_momenta.empty()) std::cout << " <none>";
    std::cout << '\n';

    std::cout << "  Rayleigh momenta =";
    for(const auto& rule : this->Rayleigh_momenta)
      {
        std::cout << " " << rule.first << " -> " << rule.second << ";";
      }
    if(this->Rayleigh_momenta.empty()) std::cout << " <none>";
    std::cout << '\n';
  }
