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


#include "Bk_tree.h"


Bk_tree::Bk_tree(const Bk_tree& obj)
  : loc(obj.loc),
    k1(obj.k1),
    k2(obj.k2),
    k3(obj.k3),
    name(obj.name),
    tag(obj.tag),
    Btree(obj.Btree)
  {
  }


Bk_tree& Bk_tree::operator+=(const Bk_tree& obj)
  {
    // check that bispectra use compatible momentum variables
    if(this->k1 != obj.k1) throw exception(ERROR_CANT_ADD_BK_INCOMPATIBLE_MOMENTA, exception_code::Bk_error);
    if(this->k2 != obj.k2) throw exception(ERROR_CANT_ADD_BK_INCOMPATIBLE_MOMENTA, exception_code::Bk_error);
    if(this->k3 != obj.k3) throw exception(ERROR_CANT_ADD_BK_INCOMPATIBLE_MOMENTA, exception_code::Bk_error);

    // merge tree databases
    this->Btree += obj.Btree;

    return *this;
  }


void Bk_tree::simplify(const GiNaC::exmap& map)
  {
    this->Btree.simplify(map);
  }


void Bk_tree::canonicalize_external_momenta()
  {
    this->Btree.canonicalize_external_momenta();
  }


void Bk_tree::write(std::ostream& out) const
  {
    out << LABEL_BK_TREE << ": '" << this->name << "'" << '\n';
    out << LABEL_BK_111 << '\n';
    out << this->Btree << '\n';
  }


void Bk_tree::write_Mathematica(std::ostream& out) const
  {
    this->Btree.write_Mathematica(out, this->tag + "Tree");
  }


std::ostream& operator<<(std::ostream& str, const Bk_tree& obj)
  {
    obj.write(str);
    return str;
  }


Bk_tree operator+(const Bk_tree& a, const Bk_tree& b)
  {
    Bk_tree c{a};
    c += b;

    c.name = a.name + std::string{"+"} + b.name;
    c.tag = a.tag + b.tag;

    return c;
  }
