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


#include "oneloop_db.h"


oneloop_db::oneloop_db(const oneloop_db& obj)
  {
    // add elements from obj
    this->operator+=(obj);
  }


void oneloop_db::write(std::ostream& out) const
  {
    size_t count = 0;

    for(const auto& item : this->db)
      {
        const oneloop_expression& lp = *item.second.first;
        const std::unique_ptr<oneloop_reduced_integral>& ri = item.second.second;

        std::cout << "Element " << count << "." << '\n';
        std::cout << lp << '\n';

        if(ri)
          {
            std::cout << *ri << '\n';
          }

        ++count;
      }

    if(count == 0) out << "<no elements>" << '\n';
  }


void oneloop_db::reduce_angular_integrals(service_locator& loc, bool symmetrize)
  {
    // walk through each subintegral in turn, performing angular reduction on it
    // the 'symmetrize' flag allows optional symmetrization of the loop and Rayleigh integrals
    // to accommodate 22-type integrations
    for(auto& item : this->db)
      {
        const oneloop_expression& lp = *item.second.first;
        std::unique_ptr<oneloop_reduced_integral>& ri = item.second.second;
        ri = std::make_unique<oneloop_reduced_integral>(lp, loc, symmetrize);    // will release any previous assignment
      }
  }


void oneloop_db::simplify(const GiNaC::exmap& map)
  {
    // walk through each expression, applying simplification map to the bare expression
    // and to the reduced integral if it exists
    for(auto& item : this->db)
      {
        std::unique_ptr<oneloop_expression>& expr = item.second.first;
        std::unique_ptr<oneloop_reduced_integral>& ri = item.second.second;

        expr->simplify(map);
        if(ri) ri->simplify(map);
      }
  }


void oneloop_db::canonicalize_external_momenta()
  {
    for(auto& item : this->db)
      {
        std::unique_ptr<oneloop_expression>& expr = item.second.first;
        std::unique_ptr<oneloop_reduced_integral>& ri = item.second.second;

        expr->canonicalize_external_momenta();
        if(ri) ri->canonicalize_external_momenta();
      }
  }


GiNaC::ex oneloop_db::get_UV_limit(unsigned int order) const
  {
    GiNaC::ex expr{0};

    for(const auto& item : this->db)
      {
        const std::unique_ptr<oneloop_reduced_integral>& ri = item.second.second;
        if(ri) expr += ri->get_UV_limit(order);
      }

    return GiNaC::collect_common_factors(expr);
  }


void oneloop_db::prune()
  {
    for(auto& record : this->db)
      {
        oneloop_pair& pair = record.second;
        const std::unique_ptr<oneloop_reduced_integral>& ri = pair.second;

        if(ri) ri->prune();
      }
  }


void oneloop_db::emplace(std::unique_ptr<oneloop_expression> elt)
  {
    oneloop_expression_key key{*elt};

    // determine whether an entry with this key already exists
    auto it = this->db.find(key);

    // if a matching element already exists then we should just sum up the kernels
    // (at this stage we assume that no angular reduction has been done; complain if it has)
    if(it != this->db.end())
      {
        oneloop_pair& record = it->second;

        std::unique_ptr<oneloop_expression>& kernel = record.first;
        std::unique_ptr<oneloop_reduced_integral>& reduced = record.second;

        if(reduced) throw exception(ERROR_ONELOOP_INTEGRAL_INSERT_AFTER_REDUCTION, exception_code::Pk_error);

        *kernel += *elt;
        return;
      }

    // otherwise, we need to insert a new element
    auto res = this->db.emplace(std::move(key), std::make_pair(std::move(elt), std::unique_ptr<oneloop_reduced_integral>{nullptr}));
    if(!res.second) throw exception(ERROR_ONELOOP_INTEGRAL_INSERT_FAILED, exception_code::Pk_error);
  }


void oneloop_db::clear_reduced_integrals()
  {
    for(auto& record : this->db)
      {
        oneloop_pair& pair = record.second;
        std::unique_ptr<oneloop_reduced_integral>& ri = pair.second;

        if(ri) ri.reset(nullptr);
      }
  }


oneloop_db& oneloop_db::operator+=(const oneloop_db& obj)
  {
    // walk through database from 'obj', and insert kernels as we go

    for(auto& item : obj.db)
      {
        const oneloop_pair& record = item.second;
        const std::unique_ptr<oneloop_expression>& kernel = record.first;

        if(kernel)
          {
            // copy kernel and emplace it in our own database
            this->emplace(std::make_unique<oneloop_expression>(*kernel));
          }
      }

    return *this;
  }


void oneloop_db::write_Mathematica(std::ostream& out, std::string symbol, bool do_dx) const
  {
    out << symbol << " = ";

    unsigned int count = 0;
    size_t chars_written = 0;

    for(auto& record : this->db)
      {
        const oneloop_pair& pair = record.second;
        const std::unique_ptr<oneloop_reduced_integral>& ri = pair.second;

        if(ri && !ri->empty())
          {
            if(count > 0) out << " + ";

            auto output = ri->to_Mathematica(do_dx);
            out << output;

            ++count;
            chars_written += output.length();
          }
      }

    if(chars_written == 0) out << "0";

    out << ";" << '\n';
  }


std::ostream& operator<<(std::ostream& str, const oneloop_db& obj)
  {
    obj.write(str);
    return str;
  }
