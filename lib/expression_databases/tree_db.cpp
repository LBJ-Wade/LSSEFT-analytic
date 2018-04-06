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

#include "tree_db.h"

#include "localizations/messages.h"
#include "shared/exceptions.h"


tree_db::tree_db(const tree_db& obj)
  {
    // add elements from obj
    this->operator+=(obj);
  }


void tree_db::emplace(std::unique_ptr<tree_expression> elt)
  {
    tree_expression_key key{*elt};

    // determine whether an entry with this key already exists
    auto it = this->db.find(key);

    // if a matching element already exists then we should just sum up the kernels
    if(it != this->db.end())
      {
        std::unique_ptr<tree_expression>& kernel = it->second;

        *kernel += *elt;
        return;
      }

    // otherwise, we need to insert a new element
    auto res = this->db.emplace(std::move(key), std::move(elt));
    if(!res.second) throw exception(ERROR_TREE_INSERT_FAILED, exception_code::Pk_error);
  }


tree_db& tree_db::operator+=(const tree_db& obj)
  {
    // walk through database from 'obj', and copy+insert kernels as we go
    for(auto& item : obj.db)
      {
        const std::unique_ptr<tree_expression>& record = item.second;

        if(record) this->emplace(std::make_unique<tree_expression>(*record));
      }

    return *this;
  }


void tree_db::write(std::ostream& out) const
  {
    size_t count = 0;

    for(const auto& item : this->db)
      {
        const std::unique_ptr<tree_expression>& tr = item.second;

        if(tr)
          {
            std::cout << "Element " << count << "." << '\n';
            std::cout << *tr << '\n';
            ++count;
          }
      }

    if(count == 0) out << "<no elements>" << '\n';
  }


void tree_db::write_Mathematica(std::ostream& out, std::string symbol) const
  {
    out << symbol << " = ";

    unsigned int count = 0;
    size_t chars_written = 0;

    for(auto& record : this->db)
      {
        const std::unique_ptr<tree_expression>& tr = record.second;

        if(tr)
          {
            if(count > 0) out << " + ";

            auto output = tr->to_Mathematica(false);
            out << output;

            ++count;
            chars_written += output.length();
          }
      }

    if(chars_written == 0) out << "0";

    out << ";" << '\n';
  }


void tree_db::simplify(const GiNaC::exmap& map)
  {
    // walk through each expression, applying simplification map
    for(auto& item : this->db)
      {
        std::unique_ptr<tree_expression>& expr = item.second;
        expr->simplify(map);
      }
  }


void tree_db::canonicalize_external_momenta()
  {
    for(auto& item : this->db)
      {
        std::unique_ptr<tree_expression>& expr = item.second;
        expr->canonicalize_external_momenta();
      }
  }


std::ostream& operator<<(std::ostream& str, const tree_db& obj)
  {
    obj.write(str);
    return str;
  }
