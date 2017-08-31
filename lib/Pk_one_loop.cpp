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
            std::cout << "Element " << count << "." << '\n';
            std::cout << *item << '\n';

            ++count;
          }
      }


    void Pk_db::reduce_angular_integrals()
      {
        // walk through each integral in turn
        for(const auto& ker : this->db)
          {
            ker->canonicalize_momenta();
            ker->dot_products_to_cos();
          }
      }

  }   // namespace Pk_one_loop_impl


void Pk_one_loop::reduce_angular_integrals()
  {
    // apply reduction algorithm to each loop container
    this->P13.reduce_angular_integrals();
    this->P22.reduce_angular_integrals();
  };


std::ostream& operator<<(std::ostream& str, const Pk_one_loop_impl::Pk_db& obj)
  {
    obj.write(str);
    return str;
  }
