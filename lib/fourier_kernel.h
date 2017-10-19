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

#include "services/service_locator.h"

#include "utilities/hash_combine.h"
#include "utilities/GiNaC_utils.h"

#include "shared/common.h"
#include "shared/exceptions.h"
#include "shared/error.h"
#include "localizations/messages.h"

#include "SPT/time_functions.h"

#include "ginac/ginac.h"


namespace fourier_kernel_impl
  {
    
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
    
    //! multiplication by arbitrary expression: expr * kernel
    kernel operator*(const GiNaC::ex a, const kernel& b);
    
    //! multiplication by arbitrary expression: kernel * expr
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
               service_locator& lc_);

        //! alternative constructor accepts just an initial_value_set and a symbol factory refernce;
        //! sets momentum kernel and time function to unity
        kernel(initial_value_set iv_, service_locator& sl_);
        
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

        //! convert to EdS approximation in which all time-dependent factors are multiples of D_lin
        void to_EdS();
        
        
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
        
        //! cache reference to service locator
        service_locator& loc;
        
        
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
        explicit key(kernel& k);
        
        //! alternative constructor accepts explicit references
        key(const time_function& tm_, const initial_value_set& iv_);
        
        //! destructor is default
        ~key() = default;

        
        // SERVICES
        
      public:

        //! hash
        size_t hash() const;

        //! check for equality
        bool is_equal(const key& obj) const;

      protected:

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


// specialize std::hash and std::equal_to to work for key_type
namespace std
  {
    
    template<>
    struct hash<fourier_kernel_impl::key>
      {
        size_t operator()(const fourier_kernel_impl::key& obj) const
          {
            return obj.hash();
          }
      };
    
    
    template<>
    struct equal_to<fourier_kernel_impl::key>
      {
        bool operator()(const fourier_kernel_impl::key& a, const fourier_kernel_impl::key& b) const
          {
            return a.is_equal(b);
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


//! compute Galileon-2 operator
template <unsigned int N>
fourier_kernel<N> Galileon2(const fourier_kernel<N>& a);

//! compute Galileon-3 operator
template <unsigned int N>
fourier_kernel<N> Galileon3(const fourier_kernel<N>& a);


//! compute the special convective term needed at third order in the bias expansion
template <unsigned int N>
fourier_kernel<N> convective_bias_term(const fourier_kernel<N>& vp, const fourier_kernel<N>& delta);


//! kernel represents an object defined by an integral kernel and early-time
//! values for each stochastic quantity such as the density constrast \delta*_k
//! the template parameter N represents the maximum order we wish to keep
template <unsigned int N>
class fourier_kernel
  {
    
    // TYPES
    
  public:

    //! pull in key_type
    using key_type = fourier_kernel_impl::key;
    
    //! pull in kernel type
    using kernel_type = fourier_kernel_impl::kernel;
    
    //! pull in kernel_db
    using kernel_db = fourier_kernel_impl::kernel_db;

    
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
    explicit fourier_kernel(service_locator& lc_)
      : loc(lc_)
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

    //! insert a symmetrized kernel into the database
    void insert_symmetric(time_function t, initial_value_set s, GiNaC::ex K, subs_list vs);

    //! insert a unsymmetrized kernel into the database
    void insert_raw(time_function t, initial_value_set s, GiNaC::ex K, subs_list vs);
    
    
    // OPERATIONS
    
  public:
    
    //! extract list of elements of fixed order as a new Fourier kernel
    fourier_kernel order(unsigned int ord) const;

    //! get size
    size_t size() const { return this->kernels.size(); }

    //! convert to EdS approximation in which all time-dependent functions are multiples of D_lin
    fourier_kernel<N> to_EdS();

    //! exchange with a second Fourier kernel object
    void swap(fourier_kernel<N>& obj) { this->kernels.swap(obj.kernels); }

    
    // SERVICES
    
  public:
    
    //! write self to a stream
    void write(std::ostream& out) const;
    
    
    // INTERNAL DATA
    
  private:
    
    //! cache reference to service locator
    service_locator& loc;
    
    //! database of kernels
    kernel_db kernels;
    
    
    friend class service_locator;
    
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

    friend fourier_kernel Galileon2<>(const fourier_kernel& a);
    friend fourier_kernel Galileon3<>(const fourier_kernel& b);

    friend fourier_kernel convective_bias_term<>(const fourier_kernel& vp, const fourier_kernel& delta);
    
  };


//! validate that an initial value set is not empty
bool validate_ivset_nonempty(const initial_value_set& s, const GiNaC::ex& K, bool silent);

//! validate that a substitution list doesn't overlap with an initial value list
void validate_subslist(const initial_value_set& s, const subs_list& vs);

//! validate that a given kernel has the correct structure (is a scalar, is a rational function of the momenta)
void validate_structure(const GiNaC::ex& K);

//! validate that a kernel and set of initial values match (no unknown momenta in kernel)
void validate_momenta(const initial_value_set& s, const subs_list& vs, const GiNaC::ex& K,
                      const GiNaC_symbol_set& params, bool silent);

//! compute normalization factor for a given time function
//! division by this factor will place the time function into a canonical form
GiNaC::ex get_normalization_factor(const time_function& tm, service_locator& sl);

//! perform symmetrization of kernel
using symmetrization_db = std::vector< GiNaC::exmap >;
symmetrization_db build_symmetrizations(const GiNaC::ex& K, const initial_value_set& s);

//! partition a GiNaC expression into factors, one for the
//! time factor (first member of pair) and one for the integrand (second member of pair)
std::pair<GiNaC::ex, GiNaC::ex> partition_factor(const GiNaC::ex& expr, service_locator& loc);


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
    // warn if initial value set is empty
    if(!validate_ivset_nonempty(s, K, silent)) return *this;

    // ensure that the substitution list (used to specify remappings for Rayleigh momenta)
    // is of the correct format
    validate_subslist(s, vs);

    // normalize the time function, redistributing factors into the kernel if needed
    auto norm = get_normalization_factor(t, this->loc);
    t /= norm;
    K *= norm;

    // simplify index structure in K if possible
    K = simplify_index(K);

    // validate that K is structurally OK (scalar, rational)
    validate_structure(K);
    
    // validate that momentum variables used in K match those listed in the stochastic terms
    validate_momenta(s, vs, K, this->loc.get_symbol_factory().get_parameters(), silent);

    if(this->loc.get_argunent_cache().get_auto_symmetrize())
      {
        this->insert_symmetric(std::move(t), std::move(s), std::move(K), std::move(vs));
        return *this;
      }

    this->insert_raw(std::move(t), std::move(s), std::move(K), std::move(vs));
    return *this;
  }


template <unsigned int N>
void fourier_kernel<N>::insert_symmetric(time_function t, initial_value_set s, GiNaC::ex K, subs_list vs)
  {
    // build list of symmetrizations for this initial value set
    auto sym_groups = build_symmetrizations(K, s);

    // get size of symmetrization set
    size_t perms = sym_groups.size();
    auto perms_N = GiNaC::numeric(perms);

    // loop over all perms, perform symmetrization, and insert in the kernel list
    for(const auto& perm : sym_groups)
      {
        // permute K
        auto perm_K = K.subs(perm) / perms_N;

        // permute substitution list
        subs_list perm_vs;
        for(const auto& v : vs)
          {
            perm_vs[v.first] = v.second.subs(perm);
          }

        this->insert_raw(t, s, std::move(perm_K), std::move(perm_vs));
      }
  }


template <unsigned int N>
void fourier_kernel<N>::insert_raw(time_function t, initial_value_set s, GiNaC::ex K, subs_list vs)
  {
    // construct a key for this record
    auto ker = std::make_unique<kernel_type>(std::move(K), std::move(s), std::move(t), std::move(vs), this->loc);
    key_type key{*ker};

    // now need to insert this kernel into the database; first, check whether an entry with this
    // key already exists
    auto it = this->kernels.find(key);

    // if so, then a matching element already exists and we should add the current kernel to it
    // notice that momentum relabelling, if needed, is handled by the kernel addition implementation
    if(it != this->kernels.end())
      {
        *it->second += *ker;
        return;
      }

    // otherwise, we can insert directly
    auto res = this->kernels.emplace(key, std::move(ker));
    if(!res.second) throw exception(ERROR_KERNEL_INSERT_FAILED, exception_code::kernel_error);
  }


template <unsigned int N>
fourier_kernel<N> fourier_kernel<N>::order(unsigned int ord) const
  {
    auto r = this->loc.template make_fourier_kernel<N>();
    
    if(ord > N) return std::move(r);
    
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
    
    return std::move(r);
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
    auto r = a.loc.template make_fourier_kernel<N>();
    
    // copy elements from a into this blank kernel, applying op as we go
    for(const auto& t : a.kernels)
      {
        const auto& old_ker = *t.second;
        
        // apply operator to kernel
        auto new_ker = op(old_ker);
        
        // insert new kernel
        r.add(new_ker, true);
      }
    
    return std::move(r);
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
    
    return std::move(r);
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
    return a * (GiNaC::numeric{1}/b);
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


//! generic algorithm to construct a product of two kernels, with a parametrizable rule
//! for constructing the product at each order
template <unsigned int N, typename InsertionRule>
void KernelProduct(const fourier_kernel<N>& a, const fourier_kernel<N>& b, InsertionRule rule)
  {
    // work through all orders that can appear in the product;
    // since everything is a perturbation, the product of two perturbations is second order
    // and we should start at 2

    // i labels order in product
    for(unsigned int i = 2; i <= N; ++i)
      {
        // j labels order to extract from first factor
        for(unsigned int j = 1; j <= i-1; ++j)
          {
            // extract terms of order j from a and order (i-j) from b
            auto a_set = a.order(j);
            auto b_set = b.order(i-j);
            
            // cross-multiply all of these terms and insert into r
            for(auto ta = a_set.cbegin(); ta != a_set.cend(); ++ta)
              {
                for(auto tb = b_set.cbegin(); tb != b_set.cend(); ++tb)
                  {
                    // apply InsertionRule to these factors
                    // the rule may customize the kernel before inserting it in a destination container
                    rule(*ta->second, *tb->second);
                  }
              }
          }
      }
  }


//! generic algorithm to construct a product of three kernels, with a parametrizable rule
//! for constructing the product at each order
template <unsigned int N, typename InsertionRule>
void KernelProduct(const fourier_kernel<N>& a, const fourier_kernel<N>& b, const fourier_kernel<N>& c,
                   InsertionRule rule)
  {
    // i labels order in product
    for(unsigned int i = 3; i <= N; ++i)
      {
        // j1 labels order to extract from first factor
        for(unsigned int j1 = 1; j1 <= i-2; ++j1)
          {
            // j2 labels order to extract from second factor
            for(unsigned int j2 = 1; j2 <= i-j1-1; ++j2)
              {
                // extract terms
                auto a_set = a.order(j1);
                auto b_set = b.order(j2);
                auto c_set = c.order(i-j1-j2);

                // cross-multiply all these terms and insert into r
                for(auto ta = a_set.cbegin(); ta != a_set.cend(); ++ta)
                  {
                    for(auto tb = b_set.cbegin(); tb != b_set.cend(); ++tb)
                      {
                        for(auto tc = c_set.cbegin(); tc != c_set.cend(); ++tc)
                          {
                            rule(*ta->second, *tb->second, *tc->second);
                          }
                      }
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
    auto r = a.loc.template make_fourier_kernel<N>();

    // insert step is trivial and should just copy each kernel product into the new container r
    auto ins = [&](const kernel& c, const kernel& d) -> void
      {
        auto ker = c*d;
        r.add(ker, true);
      };
    
    KernelProduct(a, b, ins);
    
    return std::move(r);
  }


template <unsigned int N>
fourier_kernel<N> diff_t(const fourier_kernel<N>& a)
  {
    const auto& z = a.loc.get_symbol_factory().get_z();
    
    return -FRW::Hub(z) * (GiNaC::numeric{1}+z) * diff_z(a);
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
        auto& sf = a.loc.get_symbol_factory();
        auto label_sym = sf.make_unique_Rayleigh_momentum();
        vector label = sf.make_vector(label_sym);

        auto rsq = -label.norm_square();
        b.multiply_kernel(GiNaC::ex(1)/rsq, label_sym, vsum.get_expr());

        return b;
      });
  }


template <unsigned int N>
fourier_kernel<N> gradgrad(const fourier_kernel<N>& a, const fourier_kernel<N>& b)
  {
    using fourier_kernel_impl::kernel;

    auto& sf = a.loc.get_symbol_factory();
    
    // manufacture a blank fourier kernel of max order N
    auto r = a.loc.template make_fourier_kernel<N>();
    
    // insert step should multiply each product kernel by -ka.kb
    auto ins = [&](kernel c, kernel d) -> void
      {
        auto idx = sf.make_unique_index();
        
        auto kc = c.get_total_momentum();
        auto kc_idx = kc[idx];
        
        auto kd = d.get_total_momentum();
        auto kd_idx = -kd[idx];
        
        // this implementation relies on kernels not noticing the dangling index during multiplication
        c.multiply_kernel(kc_idx);
        d.multiply_kernel(kd_idx);
        
        auto ker = c*d;
        r.add(ker, true);
      };
    
    KernelProduct(a, b, ins);
    
    return std::move(r);
  }


template <unsigned int N>
fourier_kernel<N> dotgrad(const vector& a, const fourier_kernel<N>& b)
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;

    // insert step should dot each product kernel with i a.kb
    return transform_kernel(b, [&](kernel c) -> kernel
      {
        auto kc = c.get_total_momentum();
        
        c.multiply_kernel(GiNaC::I*dot(a, kc));
        
        return c;
      });
  }


template <unsigned int N>
fourier_kernel<N> convective_bias_term(const fourier_kernel<N>& vp, const fourier_kernel<N>& delta)
  {
    static_assert(N == 3, "convective_bias_term() is currently configured only for N=3 kernels");

    using fourier_kernel_impl::kernel;

    auto& sf = delta.loc.get_symbol_factory();

    // manufacture a blank fourier kernel of max order N
    auto r = delta.loc.template make_fourier_kernel<N>();

    // get linear kernels from vp and delta
    auto vp_1 = vp.order(1);
    auto delta_1 = delta.order(1);

    for(auto ta = vp_1.cbegin(); ta != vp_1.cend(); ++ta)
      {
        for(auto tb = vp_1.cbegin(); tb != vp_1.cend(); ++tb)
          {
            for(auto tc = delta_1.cbegin(); tc != delta_1.cend(); ++tc)
              {
                auto idx_i = sf.make_unique_index();
                auto idx_j = sf.make_unique_index();

                auto ker_a = *ta->second;
                auto ker_b = *tb->second;
                auto ker_c = *tc->second;

                auto k_a = ker_a.get_total_momentum();
                auto k_b = ker_b.get_total_momentum();
                auto k_c = ker_c.get_total_momentum();

                auto k_a_i = k_a[idx_i];
                auto k_b_i = k_b[idx_i];
                auto k_b_j = k_b[idx_j];
                auto k_c_j = k_c[idx_j];

                auto ker_d = ker_a;
                auto ker_e = ker_b;
                auto ker_f = ker_c;

                auto k_d = ker_d.get_total_momentum();
                auto k_e = ker_e.get_total_momentum();
                auto k_f = ker_f.get_total_momentum();

                auto k_d_i = k_d[idx_i];
                auto k_e_j = k_e[idx_j];
                auto k_f_i = k_f[idx_i];
                auto k_f_j = k_f[idx_j];

                ker_a.multiply_kernel(GiNaC::I*k_a_i);
                ker_b.multiply_kernel(-k_b_i*k_b_j);
                ker_c.multiply_kernel(GiNaC::I*k_c_j);

                ker_d.multiply_kernel(GiNaC::I*k_d_i);
                ker_e.multiply_kernel(GiNaC::I*k_e_j);
                ker_f.multiply_kernel(-k_f_i*k_f_j);

                auto ker = ker_a*ker_b*ker_c + ker_d*ker_e*ker_f;
                r.add(ker, true);
              }
          }
      }

    return std::move(r);
  }


template <unsigned int N>
fourier_kernel<N> Galileon2(const fourier_kernel<N>& a)
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;

    auto& sf = a.loc.get_symbol_factory();

    // manufacture a blank fourier kernel of max order N
    auto r = a.loc.template make_fourier_kernel<N>();

    auto ins = [&](kernel c, kernel d) -> void
      {
        auto idx_i = sf.make_unique_index();
        auto idx_j = sf.make_unique_index();

        auto kc = c.get_total_momentum();
        auto kc_i = kc[idx_i];
        auto kc_j = kc[idx_j];

        auto kd = d.get_total_momentum();
        auto kd_i = kd[idx_i];
        auto kd_j = kd[idx_j];

        auto e = c;
        auto f = d;

        auto ke = e.get_total_momentum();
        auto kf = f.get_total_momentum();

        c.multiply_kernel(-kc_i*kc_j);
        d.multiply_kernel(-kd_i*kd_j);

        e.multiply_kernel(-ke.norm_square());
        f.multiply_kernel(-kf.norm_square());

        auto ker = c*d - e*f;
        r.add(ker, true);
      };

    KernelProduct(a, a, ins);

    return std::move(r);
  }


template <unsigned int N>
fourier_kernel<N> Galileon3(const fourier_kernel<N>& a)
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;

    auto& sf = a.loc.get_symbol_factory();

    // manufacture a blank fourier kernel of max order N
    auto r = a.loc.template make_fourier_kernel<N>();

    // insert step should dot each product kernel with i a.kb
    auto ins = [&](kernel c, kernel d, kernel e) -> void
      {
        auto idx_i = sf.make_unique_index();
        auto idx_j = sf.make_unique_index();
        auto idx_k = sf.make_unique_index();

        auto kc = c.get_total_momentum();
        auto kc_i = kc[idx_i];
        auto kc_j = kc[idx_j];

        auto kd = d.get_total_momentum();
        auto kd_j = kd[idx_j];
        auto kd_k = kd[idx_k];

        auto ke = e.get_total_momentum();
        auto ke_k = ke[idx_k];
        auto ke_i = ke[idx_i];

        auto f = c;
        auto g = d;
        auto h = e;

        auto kf = f.get_total_momentum();
        auto kg = g.get_total_momentum();
        auto kh = h.get_total_momentum();

        auto i = c;
        auto j = d;
        auto k = e;

        auto ki = i.get_total_momentum();
        auto ki_i = ki[idx_i];
        auto ki_j = ki[idx_j];

        auto kj = j.get_total_momentum();
        auto kj_i = kj[idx_i];
        auto kj_j = kj[idx_j];

        auto kk = k.get_total_momentum();

        c.multiply_kernel(-kc_i*kc_j);
        d.multiply_kernel(-kd_j*kd_k);
        e.multiply_kernel(-ke_k*ke_i);

        f.multiply_kernel(-kf.norm_square());
        g.multiply_kernel(-kg.norm_square());
        h.multiply_kernel(-kh.norm_square());

        i.multiply_kernel(-ki_i*ki_j);
        j.multiply_kernel(-kj_i*kj_j);
        k.multiply_kernel(-kk.norm_square());

        auto ker = -(2*c*d*e + f*g*h - 3*i*j*k)/2;
        r.add(ker, true);
      };

    KernelProduct(a, a, a, ins);

    return std::move(r);
  }


template <unsigned int N>
fourier_kernel<N> fourier_kernel<N>::to_EdS()
  {
    using fourier_kernel_impl::kernel;
    using fourier_kernel_impl::transform_kernel;

    // convert all kernels in the database and insert them into the new kernel
    return transform_kernel(*this, [&](kernel c) -> kernel
      {
        // convert time factor to EdS approximation
        c.to_EdS();
        return c;
      });
  }


#endif //LSSEFT_ANALYTIC_FOURIER_KERNEL_H
