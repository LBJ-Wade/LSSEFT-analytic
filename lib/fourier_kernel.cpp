//
// Created by David Seery on 21/08/2017.
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

#include <algorithm>
#include <sstream>

#include "fourier_kernel.h"


void validate_structure(const GiNaC::ex& K)
  {
    auto idcs = K.get_free_indices();
    if(!idcs.empty())
      {
        throw exception(ERROR_KERNEL_NOT_SCALAR, exception_code::kernel_error);
      }
    
    // validate that the kernel K is a rational function
    // more general kernels are not supported (yet)
    if(!is_rational(K))
      {
        std::cerr << K << '\n';
        throw exception(ERROR_KERNEL_NOT_RATIONAL, exception_code::kernel_error);
      }
  }


void validate_momenta(const initial_value_set& s, const GiNaC::ex& K, bool silent)
  {
    const auto used = get_expr_symbols(K);
    const auto& avail = s.get_momenta();
    
    GiNaC_symbol_set used_not_avail;
    GiNaC_symbol_set avail_not_used;
    
    std::set_difference(used.cbegin(), used.cend(), avail.cbegin(), avail.cend(),
                        std::inserter(used_not_avail, used_not_avail.end()), GiNaC_utils_impl::compare_symbol{});
    std::set_difference(avail.cbegin(), avail.cend(), used.cbegin(), used.cend(),
                        std::inserter(avail_not_used, avail_not_used.end()), GiNaC_utils_impl::compare_symbol{});
    
    // raise exception on attempt to use a Fourier momentum that isn't available
    if(!used_not_avail.empty())
      {
        std::ostringstream msg;
        msg << (used_not_avail.size() == 1 ? ERROR_UNKNOWN_MOMENTA_SING : ERROR_UNKNOWN_MOMENTA_PLURAL) << " ";
        
        // attach unknown momenta to error message
        unsigned int count = 0;
        for(const auto& k : used_not_avail)
          {
            if(count > 0) msg << ", ";
            msg << k;
            ++count;
          }
        
        // throw
        throw exception(msg.str(), exception_code::kernel_error);
      }
    
    // issue warning if there are available momenta that aren't used in the kernel
    if(!avail_not_used.empty() && !silent)
      {
        std::ostringstream msg;
        msg << (avail_not_used.size() == 1 ? WARNING_UNUSED_MOMENTA_SING : WARNING_UNUSED_MOMENTA_PLURAL) << " ";
        
        // attach unused momenta to message
        unsigned int count = 0;
        for(const auto& k : avail_not_used)
          {
            if(count > 0) msg << ", ";
            msg << k;
            ++count;
          }
        
        // issue error message
        error_handler err;
        err.warn(msg.str());
        
        std::ostringstream ker;
        ker << WARNING_KERNEL_EXPRESSION << " = " << K;
        err.info(ker.str());
      }
  }
