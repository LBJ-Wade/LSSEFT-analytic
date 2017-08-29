//
// Created by David Seery on 25/08/2017.
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

#include "Pk_one_loop.h"


namespace Pk_one_loop_impl
  {
    
    
    std::ostream& operator<<(std::ostream& str, const Pk_element& obj)
      {
        obj.write(str);
        return str;
      }
    
    
    std::ostream& operator<<(std::ostream& str, const Pk_db& obj)
      {
        size_t count = 0;
        for(const auto& item : obj)
          {
            std::cout << "Element " << count << "." << '\n';
            std::cout << *item << '\n';
            
            ++count;
          }
        
        return str;
      }
    
    
    void Pk_element::write(std::ostream& out) const
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
    
  }   // namespace Pk_one_loop_impl
