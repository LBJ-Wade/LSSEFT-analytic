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


void validate_momenta(const initial_value_set& s, const GiNaC::ex& K)
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
  }


fourier_kernel& fourier_kernel::add(time_function t, initial_value_set s, GiNaC::ex K)
  {
    // construct key for this kernel
    key_type key = std::make_pair(t, s);
    
    // validate that K is structurally OK (scalar, rational)
    validate_structure(K);
    
    // validate that momentum variables used in K match those listed in the stochastic terms
    validate_momenta(s, K);
    
    // now need to insert this kernel into the database; first, check whether an entry with this
    // key already exists
    auto it = this->kernels.find(key);
    
    // it not, we can insert directly
    if(it == this->kernels.end())
      {
        auto r = this->kernels.emplace(std::move(key), std::move(K));
        if(!r.second)
          {
            throw exception(ERROR_KERNEL_INSERT_FAILED, exception_code::kernel_error);
          }
        return *this;
      }
    
    // if so, we need to reallocate momenta in K to match the existing version
    // do this by pairing up symbols between lists in lexicographical order (we already know
    // that the symbol lists agree), and compare the momenta pointed to by each
    
    using iter_pair = std::pair< std::string, GiNaC::symbol >;
    using iter_set = std::vector< iter_pair >;
    
    iter_set existing_string;
    iter_set our_string;
    
    for(auto u = it->first.second.value_cbegin(); u != it->first.second.value_cend(); ++u)
      {
        existing_string.emplace_back(u->get_symbol().get_name(), u->get_momentum());
      }
    
    for(auto u = s.value_cbegin(); u != s.value_cend(); ++u)
      {
        our_string.emplace_back(u->get_symbol().get_name(), u->get_momentum());
      }
    
    std::sort(existing_string.begin(), existing_string.end(),
              [](const iter_pair& a, const iter_pair& b) -> bool { return a.first < b.first; });
    std::sort(our_string.begin(), our_string.end(),
              [](const iter_pair& a, const iter_pair& b) -> bool { return a.first < b.first; });
    
    GiNaC::exmap replace_rules;
    
    auto u = existing_string.cbegin();
    auto v = our_string.cbegin();
    for(; u != existing_string.cend() && v != our_string.cend(); ++u, ++v)
      {
        if(u->second != v->second)
          {
            replace_rules[v->second] = u->second;
          }
      }
    
    GiNaC::ex m;
    if(replace_rules.empty())
      {
        // no replacement to be done
        m = simplify_index(it->second + K);
      }
    else
      {
        m = simplify_index(it->second + K.subs(replace_rules));
      }
    it->second = m;
    
    return *this;
  }


void fourier_kernel::write(std::ostream& out) const
  {
    unsigned int count = 0;
    
    for(const auto& t : this->kernels)
      {
        const key_type& key = t.first;
        const GiNaC::ex& K = t.second;
        
        const time_function& tm = key.first;
        const initial_value_set& ivs = key.second;
        
        out << count << "." << '\n';
        out << "  time factor = " << tm << '\n';

        out << "  IC set =";
        for(auto u = ivs.value_cbegin(); u != ivs.value_cend(); ++u)
          {
            out << " " << u->get_symbol() << "(" << u->get_momentum() << ")";
          }
        out << '\n';
        
        out << "  kernel = " << K << '\n';
        
        ++count;
      }
  }
