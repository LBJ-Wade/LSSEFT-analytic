//
// Created by David Seery on 25/08/2017.
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

#include "Pk_one_loop.h"


namespace Pk_one_loop_impl
  {

    void Pk_db::write(std::ostream& out) const
      {
        size_t count = 0;

        for(const auto& item : this->db)
          {
            const loop_integral& lp = *item.second.first;
            const std::unique_ptr<one_loop_reduced_integral>& ri = item.second.second;

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


    void Pk_db::reduce_angular_integrals(service_locator& loc, bool symmetrize)
      {
        // walk through each subintegral in turn, performing angular reduction on it
        // the 'symmetrize' flag allows optional symmetrization of the loop and Rayleigh integrals
        // to accommodate 22-type integrations
        for(auto& item : this->db)
          {
            const loop_integral& lp = *item.second.first;
            std::unique_ptr<one_loop_reduced_integral>& ri = item.second.second;
            ri = std::make_unique<one_loop_reduced_integral>(lp, loc, symmetrize);    // will release any previous assignment
          }
      }


    void Pk_db::simplify(const GiNaC::exmap& map)
      {
        // walk through each subintegral, applying simplification map to reduced integral if it exists
        for(auto& item : this->db)
          {
            std::unique_ptr<one_loop_reduced_integral>& ri = item.second.second;
            if(ri) ri->simplify(map);
          }
      }


    void Pk_db::canonicalize_external_momenta()
      {
        for(auto& item : this->db)
          {
            std::unique_ptr<one_loop_reduced_integral>& ri = item.second.second;
            if(ri) ri->canonicalize_external_momenta();
          }
      }


    GiNaC::ex Pk_db::get_UV_limit(unsigned int order) const
      {
        GiNaC::ex expr{0};

        for(const auto& item : this->db)
          {
            const std::unique_ptr<one_loop_reduced_integral>& ri = item.second.second;
            if(ri) expr += ri->get_UV_limit(order);
          }

        return GiNaC::collect_common_factors(expr);
      }


    void Pk_db::prune()
      {
        for(auto& record : this->db)
          {
            loop_pair& pair = record.second;
            const std::unique_ptr<one_loop_reduced_integral>& ri = pair.second;

            if(ri) ri->prune();
          }
      }


    void Pk_db::emplace(std::unique_ptr<loop_integral> elt)
      {
        loop_integral_key key{*elt};

        // determine whether an entry with this key already exists
        auto it = this->db.find(key);

        // if a matching element already exists then we should just sum up the kernels
        if(it != this->db.end())
          {
            *it->second.first += *elt;
            return;
          }

        // otherwise, we need to insert a new element
        auto res = this->db.emplace(std::move(key), std::make_pair(std::move(elt), std::unique_ptr<one_loop_reduced_integral>{nullptr}));
        if(!res.second) throw exception(ERROR_LOOP_INTEGRAL_INSERT_FAILED, exception_code::Pk_error);
      }


    void Pk_db::write_Mathematica(std::ostream& out, std::string symbol, bool do_dx) const
      {
        out << symbol << " = ";

        unsigned int count = 0;
        size_t chars_written = 0;

        for(auto& record : this->db)
          {
            const loop_pair& pair = record.second;
            const std::unique_ptr<one_loop_reduced_integral>& ri = pair.second;

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

  }   // namespace Pk_one_loop_impl


void Pk_one_loop::simplify(const GiNaC::exmap& map)
  {
    this->Ptree.simplify(map);
    this->P13.simplify(map);
    this->P22.simplify(map);
  }


void Pk_one_loop::canonicalize_external_momenta()
  {
    this->Ptree.canonicalize_external_momenta();
    this->P13.canonicalize_external_momenta();
    this->P22.canonicalize_external_momenta();
  }


void Pk_one_loop::write(std::ostream& out) const
  {
    out << LABEL_PK << ": '" << this->name << "'" << '\n';
    out << LABEL_PK_TREE << '\n';
    out << this->get_tree() << '\n';
    out << LABEL_PK_13 << '\n';
    out << this->get_13() << '\n';
    out << LABEL_PK_22 << '\n';
    out << this->get_22() << '\n';
  }


void Pk_one_loop::write_Mathematica(std::ostream& out) const
  {
    this->Ptree.write_Mathematica(out, this->tag + "Tree", false);
    this->P13.write_Mathematica(out, this->tag + "P13", true);
    this->P22.write_Mathematica(out, this->tag + "P22", false);
  }


std::ostream& operator<<(std::ostream& str, const Pk_one_loop_impl::Pk_db& obj)
  {
    obj.write(str);
    return str;
  }


std::ostream& operator<<(std::ostream& str, const Pk_one_loop& obj)
  {
    obj.write(str);
    return str;
  }
