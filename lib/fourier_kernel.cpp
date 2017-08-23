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

#include "shared/exceptions.h"
#include "shared/error.h"
#include "localizations/messages.h"

#include "utilities/GiNaC_utils.h"


fourier_kernel& fourier_kernel::add(time_function t, initial_value_set s, GiNaC::ex K)
  {
    // construct key for this kernel
    key_type key = std::make_pair(t, s);
    
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
    
    // validate that momentum variables used in K match those listed in the stochastic terms
    auto used = get_expr_symbols(K);
    auto avail = s.get_momenta();
    
    GiNaC_symbol_set used_not_avail;
    GiNaC_symbol_set avail_not_used;

    std::set_difference(used.begin(), used.end(), avail.begin(), avail.end(),
                        std::inserter(used_not_avail, used_not_avail.end()));
    std::set_difference(avail.begin(), avail.end(), used.begin(), used.end(),
                        std::inserter(avail_not_used, avail_not_used.end()));
    
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
    if(!avail_not_used.empty())
      {
        std::ostringstream msg;
        msg << (avail_not_used.size() == 1 ? ERROR_UNUSED_MOMENTA_SING : ERROR_UNUSED_MOMENTA_PLURAL) << " ";
        
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
    
    return *this;
  }
