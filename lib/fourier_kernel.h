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
#include "utilities/GiNaC_utils.h"

#include "shared/exceptions.h"
#include "shared/error.h"
#include "localizations/messages.h"

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
            
            // to hash the initial value set, order its symbols lexicographically,
            // print them, and hash those
            std::vector<std::string> symbols;
            for(auto t = obj.second.value_cbegin(); t != obj.second.value_cend(); ++t)
              {
                symbols.push_back(t->get_symbol().get_name());
              }
            
            std::sort(symbols.begin(), symbols.end());
            
            std::ostringstream iv_string;
            for(const auto& sym : symbols)
              {
                iv_string << sym;
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
            // these are ordered in-situ lexicographically by momentum.
            // we do this by ordering their symbol names lexicographically
            // and testing for equality of those
            std::vector<std::string> a_symbols;
            std::vector<std::string> b_symbols;

            for(auto t = av.value_cbegin(); t != av.value_cend(); ++t)
              {
                a_symbols.push_back(t->get_symbol().get_name());
              }
    
            for(auto t = bv.value_cbegin(); t != bv.value_cend(); ++t)
              {
                b_symbols.push_back(t->get_symbol().get_name());
              }

            std::sort(a_symbols.begin(), a_symbols.end());
            std::sort(b_symbols.begin(), b_symbols.end());
    
            auto t = a_symbols.cbegin();
            auto u = b_symbols.cbegin();

            for(; t != a_symbols.cend() && u != b_symbols.cend(); ++t, ++u)
              {
                if(*t != *u) return false;
              }
            
            return true;
          }
      };
    
  }   // namespace std


// forward-declare class fourier_kernel
template <unsigned int N>
class fourier_kernel;

//! unary - on a Fourier kernel
template <unsigned int N>
fourier_kernel<N> operator-(const fourier_kernel<N>& a);

//! add two Fourier kernel objects
template <unsigned int N>
fourier_kernel<N> operator+(const fourier_kernel<N>& a, const fourier_kernel<N>& b);

//! subtract two Fourier kernel objects
template <unsigned int N>
fourier_kernel<N> operator-(const fourier_kernel<N>& a, const fourier_kernel<N>& b);

//! mulitply a Fourier kernel by a fixed expression, kernel * expr
template <unsigned int N>
fourier_kernel<N> operator*(const fourier_kernel<N>& a, const GiNaC::ex b);

//! multiply a Fourier kernel by a fixed expression, expr * kernel
template <unsigned int N>
fourier_kernel<N> operator*(const GiNaC::ex a, const fourier_kernel<N>& b);

//! divide a Fourier kernel by a fixed expression
template <unsigned int N>
fourier_kernel<N> operator/(const fourier_kernel<N>& a, const GiNaC::ex b);


//! kernel represents an object defined by an integral kernel and early-time
//! values for each stochastic quantity such as the density constrast \delta*_k
//! the template parameter N represents the maximum order we wish to keep
template <unsigned int N>
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
    
  protected:
    
    //! implementation: add a kernel
    fourier_kernel& add(time_function t, initial_value_set s, GiNaC::ex K, bool silent);
    
    
    // SERVICES
    
  public:
    
    //! write self to a stream
    void write(std::ostream& out) const;
    
    
    // INTERNAL DATA
    
  private:
    
    //! database of kernels
    kernel_db kernels;
    
    
    friend fourier_kernel operator-<>(const fourier_kernel& a);
    friend fourier_kernel operator+<>(const fourier_kernel& a, const fourier_kernel& b);
    friend fourier_kernel operator-<>(const fourier_kernel& a, const fourier_kernel& b);
    friend fourier_kernel operator*<>(const fourier_kernel& a, const GiNaC::ex b);
    friend fourier_kernel operator*<>(const GiNaC::ex a, const fourier_kernel& b);
    friend fourier_kernel operator/<>(const fourier_kernel& a, const GiNaC::ex b);
  
  };


void validate_structure(const GiNaC::ex& K);
void validate_momenta(const initial_value_set& s, const GiNaC::ex& K, bool silent);


template <unsigned int N>
fourier_kernel<N>& fourier_kernel<N>::add(time_function t, initial_value_set s, GiNaC::ex K)
  {
    return this->add(std::move(t), std::move(s), std::move(K), false);
  }


template <unsigned int N>
fourier_kernel<N>& fourier_kernel<N>::add(time_function t, initial_value_set s, GiNaC::ex K, bool silent)
  {
    // construct key for this kernel
    key_type key = std::make_pair(t, s);
    
    // validate that K is structurally OK (scalar, rational)
    validate_structure(K);
    
    // validate that momentum variables used in K match those listed in the stochastic terms
    validate_momenta(s, K, silent);
    
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


template <unsigned int N>
void fourier_kernel<N>::write(std::ostream& out) const
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


template <unsigned int N>
std::ostream& operator<<(std::ostream& str, const fourier_kernel<N>& obj)
  {
    obj.write(str);
    return str;
  }


template <unsigned int N>
fourier_kernel<N> operator-(const fourier_kernel<N>& a)
  {
    fourier_kernel<N> r;
    
    for(const auto& t : a.kernels)
      {
        const auto& key = t.first;
        const auto& K = t.second;
        
        const auto& tm = key.first;
        const auto& ivs = key.second;
        
        r.add(tm, ivs, -K, true);
      }
    
    return r;
  }


template <unsigned int N>
fourier_kernel<N> operator+(const fourier_kernel<N>& a, const fourier_kernel<N>& b)
  {
    fourier_kernel<N> r = a;
    
    for(const auto& t : b.kernels)
      {
        const auto& key = t.first;
        const auto& K = t.second;
        
        const auto& tm = key.first;
        const auto& ivs = key.second;
        
        r.add(tm, ivs, K, true);
      }
    
    return r;
  }


template <unsigned int N>
fourier_kernel<N> operator-(const fourier_kernel<N>& a, const fourier_kernel<N>& b)
  {
    return a + (-b);
  }


template <unsigned int N>
fourier_kernel<N> operator*(const fourier_kernel<N>& a, const GiNaC::ex b)
  {
    fourier_kernel<N> r;
    
    for(const auto& t : a.kernels)
      {
        const auto& key = t.first;
        const auto& K = t.second;
        
        const auto& tm = key.first;
        const auto& ivs = key.second;
        
        r.add(tm, ivs, b*K, true);
      }
    
    return r;
  }


template <unsigned int N>
fourier_kernel<N> operator*(const GiNaC::ex a, const fourier_kernel<N>& b)
  {
    return b*a;
  }


template <unsigned int N>
fourier_kernel<N> operator/(const fourier_kernel<N>& a, const GiNaC::ex b)
  {
    return a * (GiNaC::ex{1}/b);
  }


#endif //LSSEFT_ANALYTIC_FOURIER_KERNEL_H
