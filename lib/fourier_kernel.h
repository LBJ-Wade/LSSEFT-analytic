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

#ifndef LSSEFT_ANALYTIC_FOURIER_KERNEL_H
#define LSSEFT_ANALYTIC_FOURIER_KERNEL_H


#include <sstream>
#include <vector>
#include <set>
#include <unordered_map>

#include "initial_value.h"
#include "vector.h"

#include "utilities/hash_combine.h"

#include "ginac/ginac.h"


namespace fourier_kernel_impl
  {
    
    // set up necessary types
    
    //! type for functions of time
    using time_function = GiNaC::ex;
    
    //! type for kernel database key: kernels are stored by their time_function and initial_value_set
    //! two keys compare equal if their time functions match, and the initial value set contains the
    //! same *symbols* (doesn't have to be the same momenta)
    using key_type = std::pair< time_function, initial_value_set >;
    
    //! the kernel database is an unordered map of keys to kernel expressions
    using kernel_db = std::unordered_map< key_type, GiNaC::ex >;
    
  }   // namespace fourier_kernel_impl


// specialize std::hash and std::is_equal to work for key_type
namespace std
  {
    
    template<>
    struct hash<fourier_kernel_impl::key_type>
      {
        size_t operator()(const fourier_kernel_impl::key_type& obj) const
          {
            std::hash<std::string> string_hasher;
            
            // to hash the time expression, expand it completely and print
            std::ostringstream time_string;
            time_string << obj.first.expand();
            
            size_t h = string_hasher(time_string.str());
            
            // to hash the initial value set, print its symbol string and hash that
            std::ostringstream iv_string;
            for(auto t = obj.second.value_cbegin(); t != obj.second.value_cend(); ++t)
              {
                iv_string << t->get_symbol();
              }
            
            // combine both string hashes together
            hash_impl::hash_combine(h, iv_string.str());
            
            // return final value
            return h;
          }
      };
    
    
    template<>
    struct equal_to<fourier_kernel_impl::key_type>
      {
        bool operator()(const fourier_kernel_impl::key_type& a, const fourier_kernel_impl::key_type& b) const
          {
            const fourier_kernel_impl::time_function& at = a.first;
            const fourier_kernel_impl::time_function& bt = b.first;
            
            // test for equality of expressions
            auto rt = (at == bt);
            if(!static_cast<bool>(rt)) return false;
            
            const initial_value_set& av = a.second;
            const initial_value_set& bv = b.second;
            
            // test for equality of initial-value strings
            // note that we use only the symbols, not their momenta
            auto t = av.value_cbegin();
            auto u = bv.value_cbegin();

            for(; t != av.value_cend() && u != bv.value_cend(); ++t, ++u)
              {
                if(t->get_symbol() != u->get_symbol()) return false;
              }
            
            return true;
          }
      };
    
  }   // namespace std


//! kernel represents an object defined by an integral kernel and early-time
//! values for each stochastic quantity such as the density constrast \delta*_k
class fourier_kernel
  {
    
    // TYPES
    
  public:
    
    //! pull in time_function
    using time_function = fourier_kernel_impl::time_function;
    
  protected:
    
    //! pull in key_type
    using key_type = fourier_kernel_impl::key_type;
    
    //! pull in kernel_db
    using kernel_db = fourier_kernel_impl::kernel_db;
    
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
  
    //! constructor is default: creates an empty Fourier function
    fourier_kernel() = default;
    
    //! destructor
    ~fourier_kernel() = default;
    
    
    // KERNEL FUNCTIONS
    
  public:
    
    //! add a kernel of the form t * K(q1, q2, ..., qn) * s(q1, q2, ..., qn)
    //! where t is a time function, s is a string of stochastic initial conditions,
    //! and K is the Fourier kernel
    fourier_kernel& add(time_function t, initial_value_set s, GiNaC::ex K);
    
    
    // INTERNAL DATA
    
  private:
    
    //! database of kernels
    kernel_db kernels;
  
  };


#endif //LSSEFT_ANALYTIC_FOURIER_KERNEL_H
