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

    void Pk_db::write(std::ostream& out) const
      {
        size_t count = 0;
        for(const auto& item : this->db)
          {
            const loop_integral& lp = *item.first;

            std::cout << "Element " << count << "." << '\n';
            std::cout << lp << '\n';

            if(item.second)
              {
                const one_loop_reduced_integral& rd = *item.second;
                std::cout << rd << '\n';
              }

            ++count;
          }
      }


    void Pk_db::reduce_angular_integrals()
      {
        // walk through each subintegral in turn, performing angular reduction on it
        // walk through each integral in turn
        size_t count = 0;
        for(auto& item : this->db)
          {
            const loop_integral& lp = *item.first;
            item.second = std::make_unique<one_loop_reduced_integral>(lp);    // will release any previous assignment
          }
      }

  }   // namespace Pk_one_loop_impl


;


std::ostream& operator<<(std::ostream& str, const Pk_one_loop_impl::Pk_db& obj)
  {
    obj.write(str);
    return str;
  }
