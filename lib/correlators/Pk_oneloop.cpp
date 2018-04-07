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

#include "Pk_oneloop.h"


Pk_oneloop::Pk_oneloop(const Pk_oneloop& obj)
  : loc(obj.loc),
    k(obj.k),
    name(obj.name),
    tag(obj.tag),
    Ptree(obj.Ptree),
    P13(obj.P13),
    P22(obj.P22)
  {
    this->perform_angular_reduction();
  }


void Pk_oneloop::simplify(const GiNaC::exmap& map)
  {
    this->Ptree.simplify(map);
    this->P13.simplify(map);
    this->P22.simplify(map);
  }


void Pk_oneloop::canonicalize_external_momenta()
  {
    this->Ptree.canonicalize_external_momenta();
    this->P13.canonicalize_external_momenta();
    this->P22.canonicalize_external_momenta();
  }


void Pk_oneloop::write(std::ostream& out) const
  {
    out << LABEL_PK_ONELOOP << ": '" << this->name << "'" << '\n';
    out << LABEL_PK_11 << '\n';
    out << this->Ptree << '\n';
    out << LABEL_PK_13 << '\n';
    out << this->P13 << '\n';
    out << LABEL_PK_22 << '\n';
    out << this->P22 << '\n';
  }


void Pk_oneloop::write_Mathematica(std::ostream& out) const
  {
    this->Ptree.write_Mathematica(out, this->tag + "Tree");
    this->P13.write_Mathematica(out, this->tag + "P13", true);
    this->P22.write_Mathematica(out, this->tag + "P22", false);
  }


void Pk_oneloop::clear_angular_reductions()
  {
    this->P13.clear_reduced_integrals();
    this->P22.clear_reduced_integrals();
  }


void Pk_oneloop::perform_angular_reduction()
  {
    // perform angular reduction on integrands using Rayleigh algorithm
    this->P13.reduce_angular_integrals(loc, false);    // false = don't symmetrize q/s
    this->P22.reduce_angular_integrals(loc, true);     // true = symmetrize q/s

    // prune empty records from the reduced integrals
    this->P13.prune();
    this->P22.prune();
  }


Pk_oneloop& Pk_oneloop::operator+=(const Pk_oneloop& obj)
  {
    // check that power spectra use a compatible momentum variable
    if(this->k != obj.k) throw exception(ERROR_CANT_ADD_PK_INCOMPATIBLE_MOMENTA, exception_code::Pk_error);

    // clear angular decompositions before adding in new kernels
    this->clear_angular_reductions();

    // merge tree, P13 and P22 databases
    this->Ptree += obj.Ptree;
    this->P13 += obj.P13;
    this->P22 += obj.P22;

    // re-perform angular reductions
    this->perform_angular_reduction();

    return *this;
  }


std::ostream& operator<<(std::ostream& str, const Pk_oneloop& obj)
  {
    obj.write(str);
    return str;
  }


Pk_oneloop operator+(const Pk_oneloop& a, const Pk_oneloop& b)
  {
    Pk_oneloop c{a};
    c += b;

    c.name = a.name + std::string{"+"} + b.name;
    c.tag = a.tag + b.tag;

    return c;
  }
