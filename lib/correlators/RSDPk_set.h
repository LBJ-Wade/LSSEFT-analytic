//
// Created by David Seery on 05/04/2018.
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

#ifndef LSSEFT_ANALYTIC_RSDPK_ONELOOP_SET_H
#define LSSEFT_ANALYTIC_RSDPK_ONELOOP_SET_H


#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include <memory>

#include "Pk_oneloop.h"


using filter_list = std::vector<std::pair<GiNaC::symbol, unsigned int> >;


template <typename DatabaseType>
class RSDPk_set
  {

    // TYPES

  public:

    //! publish database type
    using db_type = DatabaseType;

  protected:

    // DatabaseType is intended to be a map (std::map or std::unordered_map) from a flyweight key class
    // to a std::unique_ptr containing the element type

    //! extract key type from map
    using key_type = typename DatabaseType::key_type;

    //! extract mapped type from map
    using mapped_type = typename DatabaseType::mapped_type;

  public:

    //! extract element type -- this is the type of the element that is stored in the database
    using element_type = typename mapped_type::element_type;


  protected:

    //! specify which mu powers to include in visitor pattern
    using visit_list = std::set<unsigned int>;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    RSDPk_set(GiNaC::symbol mu_, filter_list pt_, GiNaC_symbol_set sy_, std::string nm_, bool v);

    //! destructor
    ~RSDPk_set() = default;


    // EMPLACE TERMS

  public:

    //! accepts a one_loop_element and breaks it into individual powers of mu
    void emplace(const element_type& elt);

  protected:

    //! perform emplace on a specific database (ie. the database for mu0, mu2, ...)
    void emplace(std::unique_ptr<element_type> elt, DatabaseType& db, unsigned int mu_power);


    // TOOLS

  public:

    //! construct UV limit
    GiNaC::exvector get_UV_limit(unsigned int order = 2) const;

    //! query number of distinct time functions at each mu
    std::vector<std::vector<time_function> > get_time_functions() const;

    //! prune empty records from the database
    void prune();


    //! determine whether empty
    bool empty() const
      { return this->mu0.empty() && this->mu2.empty() && this->mu4.empty() && this->mu6.empty() && this->mu8.empty(); }


  protected:

    //! prune the database for a specific power of mu
    void prune(DatabaseType& db, unsigned int mu_power);


    // SERVICES

  public:

    //! write self to stream
    void write(std::ostream& out) const;

    //! write Mathematica script for each power of mu
    void write_Mathematica(std::ostream& out, std::string symbol, bool do_dx) const;

    //! apply visitor to the one-loop records for any given pattern of mus
    template <typename VisitorFunction>
    void visit(visit_list pattern, VisitorFunction f) const;

  private:

    //! write Mathematica script for a specific mu database
    void write_Mathematica_block(std::ostream& out, const DatabaseType& db, std::string symbol,
                                 std::string label, bool do_dx) const;


    // INTERNAL DATA

  private:

    //! group name (tree, 13, 22)
    const std::string name;

    //! cache reference to angular variable mu
    const GiNaC::symbol mu;

    //! cache filtering pattern
    const filter_list pattern;

    //! GiNaC expression representing filter pattern
    GiNaC::ex symbolic_filter;

    //! set of filtering symbols (used to kill any leftover terms)
    const GiNaC_symbol_set filter_symbols;

    //! exmap build from these symbols
    GiNaC::exmap filter_map;

    //! verbose flag
    bool verbose;


    // TERM-BY-TERM DATABASES

    //! mu^0 terms
    DatabaseType mu0;

    //! mu^2 terms
    DatabaseType mu2;

    //! mu^4 terms
    DatabaseType mu4;

    //! mu^6 terms
    DatabaseType mu6;

    //! mu^8 terms;
    DatabaseType mu8;

  };


//! stream insertion
template <typename DatabaseType>
std::ostream& operator<<(std::ostream& str, const RSDPk_set<DatabaseType>& obj);


template <typename DatabaseType>
RSDPk_set<DatabaseType>::RSDPk_set(GiNaC::symbol mu_, filter_list pt_,
                                   GiNaC_symbol_set sy_, std::string nm_, bool v)
  : name(std::move(nm_)),
    mu(std::move(mu_)),
    pattern(std::move(pt_)),
    filter_symbols(std::move(sy_)),
    symbolic_filter(1),
    verbose(v)
  {
    // set up an exmap that will set all filter_symbols to zero; we use this after pattern matching
    // to kill any remaining terms; eg. asking for the term linear in b1 from b1(1+b2) will give 1+b2, and we
    // need to set b2 to zero to get the required result '1'
    for(const auto& sym : filter_symbols)
      {
        filter_map[sym] = 0;
      }

    // build up GiNaC expression representing filter symbol string
    for(const auto& item : pattern)
      {
        symbolic_filter *= GiNaC::pow(item.first, item.second);
      }
  }


template <typename DatabaseType>
void
RSDPk_set<DatabaseType>::emplace(const element_type& elt)
  {
    auto filter = [&](const auto& e, unsigned int order) -> auto
      {
        // copy input element; assume that element_type doesn't capture any of its data by reference,
        // so what we get is a standalone copy.
        // That means we can later modify the parent container object if we require, without
        // changing what is stored here
        auto f = std::make_unique<element_type>(e);

        // filter for pattern; note that this assumes the pattern symbols are all in the integrand,
        // and none in the time function (which is certainly the intention, but things could go wrong
        // if somehow not all symbols are normalized out of the time function)
        for(const auto& pat : this->pattern)
          {
            f->filter(pat.first, pat.second);
          }

        // kill any remaining symbols (eg. if we are looking for the term linear in b1 we don't get b1 b2 contributions)
        f->simplify(this->filter_map);

        // filter for power of mu
        f->filter(this->mu, order);

        return f;
      };

    this->emplace(filter(elt, 0), this->mu0, 0);
    this->emplace(filter(elt, 2), this->mu2, 2);
    this->emplace(filter(elt, 4), this->mu4, 4);
    this->emplace(filter(elt, 6), this->mu6, 6);
    this->emplace(filter(elt, 8), this->mu8, 8);
  }


template <typename DatabaseType>
void
RSDPk_set<DatabaseType>::emplace(std::unique_ptr<element_type> elt, DatabaseType& db,
                                         unsigned int mu_power)
  {
    // nothing to do if element is empty
    if(elt->null()) return;

    if(this->verbose)
      {
        std::cout << "Storing '" << this->name << "' RSDPk_set contribution for filter pattern '"
                  << this->symbolic_filter << "' at mu^" << mu_power << '\n';
        std::cout << *elt << '\n';
      }

    typename DatabaseType::key_type key{*elt};

    // check whether a compatible entry is already present in the database
    auto it = db.find(key);

    // if such an entry is present we can just compose the integrands
    if(it != db.end())
      {
        *it->second += *elt;
        return;
      }

    // otherwise, we need to insert a new element
    auto res = db.emplace(std::move(key), std::move(elt));
    if(!res.second) throw exception(ERROR_ONELOOP_ELEMENT_INSERT_FAILED_RSD, exception_code::expression_error);
  }


template <typename DatabaseType>
GiNaC::exvector
RSDPk_set<DatabaseType>::get_UV_limit(unsigned int order) const
  {
    GiNaC::exvector values;

    auto build = [&](const DatabaseType& db) -> auto
      {
        GiNaC::ex value{0};

        for(const auto& record : db)
          {
            const auto& data = record.second;

            if(data) value += data->get_UV_limit(order);
          }

        return value;
      };

    values.push_back(build(this->mu0));
    values.push_back(build(this->mu2));
    values.push_back(build(this->mu4));
    values.push_back(build(this->mu6));
    values.push_back(build(this->mu8));

    return values;
  }


template <typename DatabaseType>
void
RSDPk_set<DatabaseType>::write(std::ostream& out) const
  {
    out << "-- mu^0" << '\n'; if(!this->mu0.empty()) out << this->mu0 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^2" << '\n'; if(!this->mu2.empty()) out << this->mu2 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^4" << '\n'; if(!this->mu4.empty()) out << this->mu4 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^6" << '\n'; if(!this->mu6.empty()) out << this->mu6 << '\n'; else out << "   <empty>" << '\n';
    out << "-- mu^8" << '\n'; if(!this->mu8.empty()) out << this->mu8 << '\n'; else out << "   <empty>" << '\n';
  }


template <typename DatabaseType>
std::vector<std::vector<time_function> >
RSDPk_set<DatabaseType>::get_time_functions() const
  {
    std::vector<std::vector<time_function> > values(5);

    auto build = [](const DatabaseType& db, std::vector<time_function>& dest)
      {
        for(const auto& item : db)
          {
            const GiNaC::ex& tm = item.first.get_time_function();

            auto it = std::find_if(dest.begin(), dest.end(),
                                   [&](const time_function& t) -> bool
                                     { return static_cast<bool>(t == tm); });

            if(it == dest.end()) dest.push_back(tm);
          }
      };

    build(this->mu0, values[0]);
    build(this->mu2, values[1]);
    build(this->mu4, values[2]);
    build(this->mu6, values[3]);
    build(this->mu8, values[4]);

    return values;
  }


template <typename DatabaseType>
void
RSDPk_set<DatabaseType>::prune()
  {
    this->prune(this->mu0, 0);
    this->prune(this->mu2, 2);
    this->prune(this->mu4, 4);
    this->prune(this->mu6, 6);
    this->prune(this->mu8, 8);
  }


template <typename DatabaseType>
void
RSDPk_set<DatabaseType>::prune(DatabaseType& db, unsigned int mu_power)
  {
    auto t = db.begin();

    while(t != db.end())
      {
        const element_type& elt = *t->second;

        // remove this element if it is null
        if(elt.null())
          {
            if(this->verbose)
              {
                std::cout << "Pruning '" << this->name << "' RSDPk_set contribution for filter pattern '"
                          << this->symbolic_filter << "' at mu^" << mu_power << '\n';
                std::cout << elt << '\n';
              }
            t = db.erase(t);
          }
        else
          {
            ++t;
          }
      }
  }


template <typename DatabaseType>
void RSDPk_set<DatabaseType>::write_Mathematica(std::ostream& out, std::string symbol, bool do_dx) const
  {
    this->write_Mathematica_block(out, this->mu0, symbol, "mu0", do_dx);
    this->write_Mathematica_block(out, this->mu2, symbol, "mu2", do_dx);
    this->write_Mathematica_block(out, this->mu4, symbol, "mu4", do_dx);
    this->write_Mathematica_block(out, this->mu6, symbol, "mu6", do_dx);
    this->write_Mathematica_block(out, this->mu8, symbol, "mu8", do_dx);
  }


template <typename DatabaseType>
void
RSDPk_set<DatabaseType>::write_Mathematica_block(std::ostream& out, const DatabaseType& db, std::string symbol,
                                                 std::string label, bool do_dx) const
  {
    out << symbol << "z" << label << " = ";

    unsigned int count = 0;
    size_t chars_written = 0;

    for(auto& record : db)
      {
        const std::unique_ptr<element_type>& elt = record.second;

        if(count > 0) out << " + ";

        auto output = elt->to_Mathematica(do_dx);
        out << output;

        ++count;
        chars_written += output.length();
      }

    if(chars_written == 0) out << "0";

    out << ";" << '\n';
  }


template <typename DatabaseType>
template <typename VisitorFunction>
void RSDPk_set<DatabaseType>::visit(visit_list pattern, VisitorFunction f) const
  {
    auto t0 = pattern.find(0);
    auto t2 = pattern.find(2);
    auto t4 = pattern.find(4);
    auto t6 = pattern.find(6);
    auto t8 = pattern.find(8);

    auto visit = [&](const DatabaseType& db) -> void
      {
        for(const auto& item : db)
          {
            const element_type& elt = *item.second;
            f(elt);
          }
      };

    if(t0 != pattern.end()) visit(this->mu0);
    if(t2 != pattern.end()) visit(this->mu2);
    if(t4 != pattern.end()) visit(this->mu4);
    if(t6 != pattern.end()) visit(this->mu6);
    if(t8 != pattern.end()) visit(this->mu8);
  }


template <typename DatabaseType>
std::ostream& operator<<(std::ostream& str, const RSDPk_set<DatabaseType>& obj)
  {
    obj.write(str);
    return str;
  }


#endif //LSSEFT_ANALYTIC_RSDPK_ONELOOP_SET_H
