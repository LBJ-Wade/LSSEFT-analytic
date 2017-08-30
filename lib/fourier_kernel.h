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

#include "services/symbol_factory.h"

#include "utilities/hash_combine.h"
#include "utilities/GiNaC_utils.h"

#include "shared/exceptions.h"
#include "shared/error.h"
#include "localizations/messages.h"

#include "SPT/time_functions.h"

#include "ginac/ginac.h"


namespace fourier_kernel_impl
  {
    
    // set up necessary types

    //! list of GiNaC non-rotationally-invariant vectors that appear in the
    //! denominator of the kernel, and which will need to be computed
    //! by Rayleigh expansion in correlation functions
    using subs_list = GiNaC::exmap;
    
    //! type for functions of time
    using time_function = GiNaC::ex;
    
    
    //! forward-declare key
    class key;
    
    //! forward-declare kernel
    class kernel;
    
    //! perform stream insertion
    std::ostream& operator<<(std::ostream& str, const kernel& a);
    
    //! unary - for kernel
    kernel operator-(const kernel& a);
    
    //! addition of kernels
    kernel operator+(const kernel& a, const kernel& b);
    
    //! subtraction of kernels
    kernel operator-(const kernel& a, const kernel& b);
    
    //! multiplication by arbitrary expression
    kernel operator*(const GiNaC::ex a, const kernel& b);
    
    //! multiplication by arbitrary expression, other way ground
    kernel operator*(const kernel& a, const GiNaC::ex b);
    
    //! multiplication of kernels
    kernel operator*(const kernel& a, const kernel& b);
    
    //! division by arbitrary expression
    kernel operator/(const kernel& a, const GiNaC::ex b);
    
    //! differentiation with respect to redshift z
    kernel diff_z(const kernel& a);
    
    
    //! kernel captures the details of a Fourier kernel
    class kernel
      {
        
        // TYPES
        
      protected:
        
        //! a momenta_list is a list of initial_value references, ordered by symbol name
        using momenta_list = std::vector< std::reference_wrapper<const initial_value> >;
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! constructor accepts a GiNaC expression and a list corresponding to
        //! non-rotationally-invariant combinations of momenta that appear in the denominator;
        //! these will need to be handled by a Rayleigh plane wave expansion when computing
        //! correlation functions
        kernel(GiNaC::ex K_, initial_value_set iv_, time_function tm_, subs_list vs_,
               symbol_factory& sf_)
          : K(std::move(K_)),
            tm(std::move(tm_)),
            iv(std::move(iv_)),
            vs(std::move(vs_)),
            sf(sf_)
          {
          }

        //! alternative constructor accepts just an initial_value_set and a symbol factory refernce;
        //! sets momentum kernel and time function to unity
        kernel(initial_value_set iv_, symbol_factory& sf_)
          : K(1),
            tm(1),
            iv(std::move(iv_)),
            sf(sf_)
          {
          }
        
        //! destructor us default
        ~kernel() = default;
        
        
        // INTERFACE
        
      public:
        
        //! get order of this term
        size_t order() const { return this->iv.size(); }
        
        //! get total momentum = sum of dummy momenta
        vector get_total_momentum() const;
        
        //! allow implicit conversion to a GiNaC::ex
        explicit operator GiNaC::ex() const { return this->K; }
    
    
        // ACCESSORS
  
      public:
    
        //! get time function
        const time_function& get_time_function() const { return this->tm; }
    
        //! get initial value list
        const initial_value_set& get_initial_value_set() const { return this->iv; }
        
        //! get substitution list
        const subs_list& get_substitution_list() const { return this->vs; }
        
        //! get kernel expression
        const GiNaC::ex& get_kernel() const { return this->K; }
        
        
        // OPERATIONS
        
      public:
        
        //! allow in-place kernel addition
        kernel& operator+=(const kernel& rhs);
        
        //! allow in-place multiplication
        kernel& operator*=(const kernel& rhs);
        
        //! adjust kernel by a factor
        kernel& multiply_kernel(GiNaC::ex f);
        
        //! adjust kernel by a factor, generating a new substitution rule
        kernel& multiply_kernel(GiNaC::ex f, GiNaC::symbol s, GiNaC::ex rule);
        
        
        // SERVICES
        
      public:
        
        //! write self to stream
        void write(std::ostream& out) const;
        
        
        // INTERNAL API
        
      protected:
        
        //! get list of momentum labels, in order corresponding to lexical order of symbols
        momenta_list get_ordered_momenta() const;
        
        
        // INTERNAL DATA
        
      private:
        
        // AGENTS
        
        //! cache reference to symbol factory
        symbol_factory& sf;
        
        
        // KERNEL DATA
        
        //! cache kernel expression
        GiNaC::ex K;
        
        //! time function
        time_function tm;
        
        //! cache initial value set
        initial_value_set iv;
        
        //! cache list of non-invariant denominator combinations
        subs_list vs;
        
        
        friend class key;
        
        friend kernel operator-(const kernel& a);
        friend kernel operator*(const GiNaC::ex a, const kernel& b);
        friend kernel operator*(const kernel& a, const kernel& b);
        
        friend kernel diff_z(const kernel& a);
        
      };
    
    
    //! key is a flyweight that captures the time and initial value combinations
    //! from a given kernel, which we use to index the database
    class key
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! constructor accepts a kernel object and captures its data
        explicit key(kernel& k)
          : tm(k.tm),
            iv(k.iv)
          {
          }
        
        //! alternative constructor accepts explicit references
        key(const time_function& tm_, const initial_value_set& iv_)
          : tm(tm_),
            iv(iv_)
          {
          }
        
        //! destructor is default
        ~key() = default;
        
        
        // ACCESSORS
        
      public:
        
        //! get time function
        const time_function& get_time_function() const { return this->tm; }
        
        //! get initial value list
        const initial_value_set& get_initial_value_set() const { return this->iv; }
        
        
        // SERVICES
        
      public:
        
        //! get lexicographically-ordered list of symbols in the initial value set;
        //! used for ordering and comparison
        std::vector<GiNaC::symbol> get_ordered_iv_symbols() const;
        
        
        // INTERNAL DATA
        
      public:
        
        //! reference to time function
        const time_function& tm;
        
        //! reference to initial value set
        const initial_value_set& iv;
        
      };
    
    //! the kernel database is an unordered map of keys to kernel expressions
    using kernel_db = std::unordered_map< key, std::unique_ptr<kernel> >;
    
  }   // namespace fourier_kernel_impl


// specialize std::hash and std::is_equal to work for key_type
namespace std
  {
    
    template<>
    struct hash<fourier_kernel_impl::key>
      {
        size_t operator()(const fourier_kernel_impl::key& obj) const
          {
            std::hash<std::string> string_hasher;
            
            // to hash the time expression, expand it completely and print
            std::ostringstream time_string;
            time_string << obj.get_time_function().expand();
            
            size_t h = string_hasher(time_string.str());
            
            // to hash the initial value set, order its symbols lexicographically
            const auto symbols = obj.get_ordered_iv_symbols();
    
            std::string symbol_string;
            std::for_each(symbols.begin(), symbols.end(),
                          [&](const GiNaC::symbol& e) -> std::string
                            { return symbol_string += e.get_name(); });
            
            // and hash that. Then we combine both hashes together.
            hash_impl::hash_combine(h, symbol_string);
            
            // return final value
            return h;
          }
      };
    
    
    template<>
    struct equal_to<fourier_kernel_impl::key>
      {
        bool operator()(const fourier_kernel_impl::key& a, const fourier_kernel_impl::key& b) const
          {
            const auto& at = a.get_time_function();
            const auto& bt = b.get_time_function();
            
            // test for equality of expressions
            auto rt = (at == bt);
            if(!static_cast<bool>(rt)) return false;
            
            // test for equality of initial-value strings
            // we do this by ordering their symbol names lexicographically
            // and testing for equality of those
            auto a_symbols = a.get_ordered_iv_symbols();
            auto b_symbols = b.get_ordered_iv_symbols();
    
            return std::equal(a_symbols.cbegin(), a_symbols.cend(),
                              b_symbols.cbegin(), b_symbols.cend(),
                              [](const GiNaC::symbol& asym, const GiNaC::symbol& bsym) -> bool
                                { return asym.get_name() == bsym.get_name(); });
          }
      };
    
  }   // namespace std


// pull in 'kernel' concept since it has wider utility
using fourier_kernel_impl::kernel;


// forward-declare class fourier_kernel
template <unsigned int N>
class fourier_kernel;

//! generic kernel transformation function
namespace fourier_kernel_impl
  {

    //! build a new kernel by applying Op to the elements of kernel a
    template <unsigned int N, typename Operator>
    fourier_kernel<N> transform_kernel(const fourier_kernel<N>& a, Operator op);

    //! build a new kernel by applying Op to the elements of kernels a and b,
    //! and then inserting elements from b into a
    template <unsigned int N, typename Operator>
    fourier_kernel<N> transform_kernel(const fourier_kernel<N>& a, const fourier_kernel<N>& b, Operator op);
    
  }   // namespace fourier_kernel_impl

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

//! multiply a Fourier kernel by a fixed expression, kernel * expr
template <unsigned int N>
fourier_kernel<N> operator*(const fourier_kernel<N>& a, const GiNaC::ex b);

//! multiply two Fourier kernels
template <unsigned int N>
fourier_kernel<N> operator*(const fourier_kernel<N>& a, const fourier_kernel<N>& b);

//! compute t-derivative of a Fourier kernel
template <unsigned int N>
fourier_kernel<N> diff_t(const fourier_kernel<N>& a);

//! compute z-derivative of a Fourier kernel
template <unsigned int N>
fourier_kernel<N> diff_z(const fourier_kernel<N>& a);

//! compute Laplacian of a Fourier kernel
template <unsigned int N>
fourier_kernel<N> Laplacian(const fourier_kernel<N>& a);

//! compute inverse Laplacian of a Fourier kernel
template <unsigned int N>
fourier_kernel<N> InverseLaplacian(const fourier_kernel<N>& a);

//! compute (grad a).(grad b) for two Fourier kernels
template <unsigned int N>
fourier_kernel<N> gradgrad(const fourier_kernel<N>& a, const fourier_kernel<N>& b);

//! compute vector.(grad b) for a vector and a Fourier kernel.
//! The Fourier convention is that f(x) = \int d^3 k (2pi)^-3 f(k) exp(ik.x),
//! so the gradient will pull down a factor of +ik
template <unsigned int N>
fourier_kernel<N> dotgrad(const vector& a, const fourier_kernel<N>& b);


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
    
    //! pull in key_type
    using key_type = fourier_kernel_impl::key;
    
    //! pull in kernel type
    using kernel_type = fourier_kernel_impl::kernel;
    
    //! pull in kernel_db
    using kernel_db = fourier_kernel_impl::kernel_db;
    
    //! pull in subs_list
    using subs_list = fourier_kernel_impl::subs_list;
    
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
  
    //! default constructor is disabled
    fourier_kernel() = delete;
    
    //! destructor
    ~fourier_kernel() = default;
    
    //! copy constructor is disabled
    fourier_kernel(const fourier_kernel& obj) = delete;
    
    //! force compiler to generate default move constructor
    fourier_kernel(fourier_kernel&& obj) = default;
    
  protected:
    
    //! constructor captures symbol_factory class and creates an empty Fourier representation
    explicit fourier_kernel(symbol_factory& sf_)
      : sf(sf_)
      {
      }
    
    
    // ITERATORS
    
  public:
    
    kernel_db::const_iterator cbegin() const { return this->kernels.cbegin(); }
    kernel_db::const_iterator cend() const   { return this->kernels.cend(); }
    
    
    // KERNEL FUNCTIONS
    
  public:
    
    //! add a kernel of the form t * K(q1, q2, ..., qn) * s(q1, q2, ..., qn)
    //! where t is a time function, s is a string of stochastic initial conditions,
    //! and K is the Fourier kernel. The optional list of GiNaC expressions
    //! noninv_list can be used to indicate rotationally-noninvariant combinations that
    //! appear in the denominator, and which should be handled by Rayleigh plane-wave
    //! expansion when computing correlation functions.
    //! WARNING: if any rotationally noninvariant expression appear explicitly in the
    //! denominator then any computed correlation functions are unlikely to be reliable
    fourier_kernel&
    add(time_function t, initial_value_set s, GiNaC::ex K, subs_list vs = subs_list{});
    
    //! add a kernel
    fourier_kernel& add(kernel_type ker);
    
  protected:
    
    //! implementation: add a kernel
    fourier_kernel&
    add(time_function t, initial_value_set s, GiNaC::ex K, subs_list vs, bool silent);
    
    //! implementation: add a kernel
    fourier_kernel& add(kernel_type ker, bool silent);
    
    
    // OPERATIONS
    
  public:
    
    //! extract list of elements of fixed order as a new Fourier kernel
    fourier_kernel order(unsigned int ord) const;

    //! get size
    size_t size() const { return this->kernels.size(); }
    
    
    // SERVICES
    
  public:
    
    //! write self to a stream
    void write(std::ostream& out) const;
    
    
    // INTERNAL DATA
    
  private:
    
    //! cache reference to symbol factory
    symbol_factory& sf;
    
    //! database of kernels
    kernel_db kernels;
    
    
    friend class symbol_factory;
    
    // friend declarations can't be partial specializations, so we must use the fully templated versions here
    template <unsigned int M, typename Operator>
    friend fourier_kernel<M> fourier_kernel_impl::transform_kernel(const fourier_kernel<M>& a, Operator op);
    template <unsigned int M, typename Operator>
    friend fourier_kernel<M> fourier_kernel_impl::transform_kernel(const fourier_kernel<M>& a, const fourier_kernel<M>& b, Operator op);
    
    friend fourier_kernel operator-<>(const fourier_kernel& a);
    friend fourier_kernel operator+<>(const fourier_kernel& a, const fourier_kernel& b);
    friend fourier_kernel operator-<>(const fourier_kernel& a, const fourier_kernel& b);
    friend fourier_kernel operator*<>(const fourier_kernel& a, const GiNaC::ex b);
    friend fourier_kernel operator*<>(const GiNaC::ex a, const fourier_kernel& b);
    friend fourier_kernel operator/<>(const fourier_kernel& a, const GiNaC::ex b);
    friend fourier_kernel operator*<>(const fourier_kernel& a, const fourier_kernel& b);
    
    friend fourier_kernel diff_t<>(const fourier_kernel& a);
    friend fourier_kernel diff_z<>(const fourier_kernel& a);
    
    friend fourier_kernel Laplacian<>(const fourier_kernel& a);
    friend fourier_kernel InverseLaplacian<>(const fourier_kernel& a);
    friend fourier_kernel gradgrad<>(const fourier_kernel& a, const fourier_kernel& b);
    friend fourier_kernel dotgrad<>(const vector& a, const fourier_kernel& b);
    
  };


//! validate that an initial value set is not empty
bool validate_ivset_nonempty(const initial_value_set& s, const GiNaC::ex& K, bool silent);

//! validate that a substitution list doesn't overlap with an initial value list
void validate_subslist(const initial_value_set& s, const fourier_kernel_impl::subs_list& vs);

//! validate that a given kernel has the correct structure (is a scalar, is a rational function of the momenta)
void validate_structure(const GiNaC::ex& K);

//! validate that a kernel and set of initial values match (no unknown momenta in kernel)
void
validate_momenta(const initial_value_set& s, const fourier_kernel_impl::subs_list& vs, const GiNaC::ex& K,
                 const GiNaC_symbol_set& params, bool silent);


template <unsigned int N>
fourier_kernel<N>& fourier_kernel<N>::add(kernel_type k)
  {
    return this->add(k.get_time_function(), k.get_initial_value_set(), k.get_kernel(), k.get_substitution_list(), false);
  }


template <unsigned int N>
fourier_kernel<N>& fourier_kernel<N>::add(kernel_type k, bool silent)
  {
    return this->add(k.get_time_function(), k.get_initial_value_set(), k.get_kernel(), k.get_substitution_list(), silent);
  }


template <unsigned int N>
fourier_kernel<N>&
fourier_kernel<N>::add(time_function t, initial_value_set s, GiNaC::ex K, subs_list vs)
  {
    return this->add(std::move(t), std::move(s), std::move(K), std::move(vs), false);
  }


template <unsigned int N>
fourier_kernel<N>&
fourier_kernel<N>::add(time_function t, initial_value_set s, GiNaC::ex K, subs_list vs, bool silent)
  {
    if(!validate_ivset_nonempty(s, K, silent)) return *this;
    validate_subslist(s, vs);
    
    // validate that K is structurally OK (scalar, rational)
    validate_structure(K);
    
    // validate that momentum variables used in K match those listed in the stochastic terms
    validate_momenta(s, vs, K, this->sf.get_parameters(), silent);
    
    // check whether an entry with this key already exists
    auto ker = std::make_unique<kernel_type>(std::move(K), std::move(s), std::move(t), std::move(vs), this->sf);
    key_type key{*ker};
    
    // now need to insert this kernel into the database; first, check whether an entry with this
    // key already exists
    auto it = this->kernels.find(key);
    
    // if so, then a matching element already exists and we should add the current kernel to it
    // notice that momentum relabelling, if needed, is handled by the kernel addition implementation
    if(it != this->kernels.end())
      {
        *it->second += *ker;
        return *this;
      }
    
    // otherwise, we can insert directly
    auto res = this->kernels.emplace(std::move(key), std::move(ker));
    if(!res.second) throw exception(ERROR_KERNEL_INSERT_FAILED, exception_code::kernel_error);

    return *this;
  }


template <unsigned int N>
fourier_kernel<N> fourier_kernel<N>::order(unsigned int ord) const
  {
    auto r = this->sf.template make_fourier_kernel<N>();
    
    if(ord > N) return r;
    
    for(const auto& t : this->kernels)
      {
        // if order matches required value, make a copy of the kernel and emplace it
        if(t.second->order() == ord)
          {
            auto copy_ker = std::make_unique<kernel_type>(*t.second);
            key_type copy_key{*copy_ker};
            auto res = r.kernels.emplace(std::move(copy_key), std::move(copy_ker));
            if(!res.second) throw exception(ERROR_KERNEL_COPY_INSERT_FAILED, exception_code::kernel_error);
          }
      }
    
    return r;
  }


template <unsigned int N>
void fourier_kernel<N>::write(std::ostream& out) const
  {
    unsigned int count = 0;
    
    for(const auto& t : this->kernels)
      {
        const auto& key = t.first;
        const auto& ker = *t.second;
        
        const time_function& tm = ker.get_time_function();
        const initial_value_set& ivs = ker.get_initial_value_set();
        const subs_list& vs = ker.get_substitution_list();
        
        out << "Kernel " << count << "." << '\n';
        t.second->write(out);
        
        ++count;
      }
  }


template <unsigned int N>
std::ostream& operator<<(std::ostream& str, const fourier_kernel<N>& obj)
  {
    obj.write(str);
    return str;
  }


template <unsigned int N, typename Operator>
fourier_kernel<N> fourier_kernel_impl::transform_kernel(const fourier_kernel<N>& a, Operator op)
  {
    // manufacture a blank fourier kernel of max order N
    auto r = a.sf.template make_fourier_kernel<N>();
    
    // copy elements from a into this blank kernel, applying op as we go
    for(const auto& t : a.kernels)
      {
        const auto& old_ker = *t.second;
        
        // apply operator to kernel
        auto new_ker = op(old_ker);
        
        // insert new kernel
        r.add(new_ker, true);
      }
    
    return r;
  }


template <unsigned int N, typename Operator>
fourier_kernel<N> fourier_kernel_impl::transform_kernel(const fourier_kernel<N>& a, const fourier_kernel<N>& b, Operator op)
  {
    // copy elements of a into a new kernel, applying op
    auto r = fourier_kernel_impl::transform_kernel(a, op);
    
    // now copy elements from b into this kernel, applying op
    for(const auto& t : b.kernels)
      {
        const auto& old_ker = *t.second;
        
        // apply operator to kernel
        auto new_ker = op(old_ker);
        
        // insert this new kernel
        r.add(new_ker, true);
      }
    
    return r;
  }


template <unsigned int N>
fourier_kernel<N> operator-(const fourier_kernel<N>& a)
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;
    
    return transform_kernel(a, [](const kernel& b) -> kernel { return -b; });
  }


template <unsigned int N>
fourier_kernel<N> operator*(const fourier_kernel<N>& a, const GiNaC::ex b)
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;
    
    return transform_kernel(a, [&](const kernel& c) -> kernel { return b*c; });
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


template <unsigned int N>
fourier_kernel<N> operator+(const fourier_kernel<N>& a, const fourier_kernel<N>& b)
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;
    
    return transform_kernel(a, b, [](const kernel& c) -> kernel { return c; });
  }


template <unsigned int N>
fourier_kernel<N> operator-(const fourier_kernel<N>& a, const fourier_kernel<N>& b)
  {
    return a + (-b);
  }


inline GiNaC::ex BasicKernelProduct(const GiNaC::ex& a, const GiNaC::ex& b,
                                    const initial_value_set& a_ivs, const initial_value_set& b_ivs)
  {
    return a*b;
  }


//! generic algorithm to construct a product of kernels, with a parametrizable rule
//! for constructing the product at a fixed order
template <unsigned int N, typename ProductRule>
void KernelProduct(const fourier_kernel<N>& a, const fourier_kernel<N>& b, ProductRule rule)
  {
    // work through all orders that can appear in the product;
    // since everything is a perturbation, the product of two perturbations is second order
    // and we should start at 2
    for(unsigned int i = 2; i <= N; ++i)
      {
        for(unsigned int j = 1; j < i; ++j)
          {
            // extract terms of order j from a and order (i-j) from b
            auto a_set = a.order(j);
            auto b_set = b.order(i-j);
            
            // cross-multiply all of these terms and insert into r
            for(auto ta = a_set.cbegin(); ta != a_set.cend(); ++ta)
              {
                for(auto tb = b_set.cbegin(); tb != b_set.cend(); ++tb)
                  {
                    // apply ProductRule to these factors
                    // the rule may customize the kernel before inserting it in a destination container
                    rule(*ta->second, *tb->second);
                  }
              }
          }
      }
  }


template <unsigned int N>
fourier_kernel<N> operator*(const fourier_kernel<N>& a, const fourier_kernel<N>& b)
  {
    using fourier_kernel_impl::kernel;
    
    // manufacture a blank fourier kernel of max order N
    auto r = a.sf.template make_fourier_kernel<N>();

    // insert step is trivial and should just copy each kernel product into the new container r
    auto ins = [&](const kernel& c, const kernel& d) -> void
      {
        auto k = c*d;
        
        r.add(k, true);
      };
    
    KernelProduct(a, b, ins);
    
    return r;
  }


template <unsigned int N>
fourier_kernel<N> diff_t(const fourier_kernel<N>& a)
  {
    const auto& z = a.sf.get_z();
    
    return -FRW::Hub(z) * (1+z) * diff_z(a);
  }


template <unsigned int N>
fourier_kernel<N> diff_z(const fourier_kernel<N>& a)
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;
    
    return transform_kernel(a, [](const kernel& b) -> kernel { return diff_z(b); });
  }


template <unsigned int N>
fourier_kernel<N> Laplacian(const fourier_kernel<N>& a)
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;
    
    return transform_kernel(a, [](kernel b) -> kernel
      {
        // extract -k^2 for this kernel
        // it's preferable to avoid introducing a new substitution rule here; those are best kept
        // for factors in the denominator
        auto vsum = b.get_total_momentum();
        auto vsq = -vsum.norm_square();
        
        b.multiply_kernel(vsq);
        return b;
      });
  }


template <unsigned int N>
fourier_kernel<N> InverseLaplacian(const fourier_kernel<N>& a)
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;
    
    return transform_kernel(a, [&](kernel b) -> kernel
      {
        // extract -k^2 for this kernel
        auto vsum = b.get_total_momentum();
        
        // generate a new substitution rule because we are introducing
        // a potentially non-rotationally invariant denominator
        auto label_sym = a.sf.make_unique_Rayleigh_momentum();
        vector label = a.sf.make_vector(label_sym);

        auto rsq = -label.norm_square();
        b.multiply_kernel(GiNaC::ex(1)/rsq, label_sym, vsum.get_expr());

        return b;
      });
  }


template <unsigned int N>
fourier_kernel<N> gradgrad(const fourier_kernel<N>& a, const fourier_kernel<N>& b)
  {
    using fourier_kernel_impl::kernel;
    
    // manufacture a blank fourier kernel of max order N
    auto r = a.sf.template make_fourier_kernel<N>();
    
    // insert step should multiply each product kernel by -ka.kb
    auto ins = [&](kernel c, kernel d) -> void
      {
        auto idx = a.sf.make_unique_index();
        
        auto kc = c.get_total_momentum();
        auto kc_idx = GiNaC::indexed(kc.get_expr(), idx);
        
        auto kd = d.get_total_momentum();
        auto kd_idx = -GiNaC::indexed(kd.get_expr(), idx);
        
        // this implementation relies on kernels not noticing the dangling index during multiplication
        c.multiply_kernel(kc_idx);
        d.multiply_kernel(kd_idx);
        
        auto k = c*d;
        
        r.add(k, true);
      };
    
    KernelProduct(a, b, ins);
    
    return r;
  }


template <unsigned int N>
fourier_kernel<N> dotgrad(const vector& a, const fourier_kernel<N>& b)
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;
    
    // manufacture a blank fourier kernel of max order N
    auto r = b.sf.template make_fourier_kernel<N>();
    
    // insert step should dot each product kernel with i a.kb
    return transform_kernel(b, [&](kernel c) -> kernel
      {
        auto kc = c.get_total_momentum();
        
        c.multiply_kernel(GiNaC::I*dot(a, kc));
        
        return c;
      });
  }



#endif //LSSEFT_ANALYTIC_FOURIER_KERNEL_H
