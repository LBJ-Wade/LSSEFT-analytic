//
// Created by David Seery on 02/10/2017.
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
#ifndef LSSEFT_ANALYTIC_LSSEFT_H
#define LSSEFT_ANALYTIC_LSSEFT_H


#include <map>
#include <functional>

#include "shared/defaults.h"

#include "lib/Pk_rsd.h"

#include "boost/filesystem/operations.hpp"


namespace LSSEFT_impl
  {

    //! holds the details of a single LSSEFT kernel
    class LSSEFT_kernel
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor
        LSSEFT_kernel(GiNaC::ex ig_, GiNaC::ex ms_, GiNaC::ex wp_, GiNaC_symbol_set vs_, GiNaC_symbol_set em_);

        //! destructor is default
        ~LSSEFT_kernel() = default;


        // ACCESSORS

      public:

        //! get external momenta
        const GiNaC_symbol_set& get_external_momenta() const { return this->external_momenta; }


        // SERVICES

      public:

        //! test for equality
        bool is_equal(const LSSEFT_kernel& obj) const;

        //! hash
        size_t hash() const;


        // FORMATTING

      public:

        //! print integrand
        std::string print_integrand(const GiNaC::exmap& subs_map) const;

        //! print measure
        std::string print_measure(const GiNaC::exmap& subs_map) const;

        //! print Wick product
        std::string print_WickProduct(const GiNaC::exmap& subs_map, const GiNaC_symbol_set& external_momenta) const;


        // INTERNAL DATA

      private:

        //! integrand
        GiNaC::ex integrand;

        //! measure
        GiNaC::ex measure;

        //! Wick product
        GiNaC::ex WickProduct;

        //! set of integration variables
        GiNaC_symbol_set variables;

        //! set of external momenta
        GiNaC_symbol_set external_momenta;

      };

  }   // namespace LSSEFT_impl


// specialize std::hash<> and std::is_equal<> to LSSEFT_impl
namespace std
  {

    template <>
    struct hash<LSSEFT_impl::LSSEFT_kernel>
      {
        size_t operator()(const LSSEFT_impl::LSSEFT_kernel& obj) const
          {
            return obj.hash();
          }
      };

    template <>
    struct equal_to<LSSEFT_impl::LSSEFT_kernel>
      {
        bool operator()(const LSSEFT_impl::LSSEFT_kernel& a, const LSSEFT_impl::LSSEFT_kernel& b) const
          {
            return a.is_equal(b);
          }
      };

  }   // namespace std


//! LSSEFT is a backend class capable of writing out C++ to implement
//! a given correlation function in the Sussex LSSEFT platform
class LSSEFT
  {

    // TYPES

  protected:

    //! power spectrum database
    using Pk_db_type = std::map< std::string, std::reference_wrapper<const Pk_rsd> >;

    //! kernel database
    using kernel_db_type = std::unordered_map< LSSEFT_impl::LSSEFT_kernel, std::string >;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    LSSEFT(boost::filesystem::path rt_, service_locator& lc_);

    //! destructor
    ~LSSEFT() = default;


    // INTERFACE

  public:

    //! add a new power spectrum
    LSSEFT& add(const Pk_rsd& P, std::string name);

    //! write output files
    void write() const;

  protected:

    //! process the kernels associated with an added power spectrum
    void process_kernels(const Pk_rsd_group& group);

    //! generate a unique kernel name
    std::string make_unique_kernel_name();


    // INTERNAL API

  protected:

    //! construct an output file name from the cached root
    boost::filesystem::path make_output_path(std::string append) const;

    //! write create block
    void write_create() const;

    //! write kernels
    void write_kernels() const;


    // INTERNAL DATA

  private:

    // SERVICES

    //! cache reference to service locator
    service_locator& loc;


    // CONFIGURATION DATA

    //! output root
    boost::filesystem::path root;

    //! current kernel number
    unsigned int kernel_count{0};

    //! kernel root string
    std::string kernel_root{LSSEFT_DEFAULT_KERNEL_ROOT};


    // DATABASES

    //! power spectrum database
    Pk_db_type Pk_db;

    //! kernel database
    kernel_db_type kernel_db;

  };


#endif //LSSEFT_ANALYTIC_LSSEFT_H
