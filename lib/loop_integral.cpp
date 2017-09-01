//
// Created by David Seery on 30/08/2017.
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

#include "loop_integral.h"

#include "detail/contractions.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


std::ostream& operator<<(std::ostream& str, const loop_integral& obj)
  {
    obj.write(str);
    return str;
  }


loop_integral::loop_integral(time_function tm_, GiNaC::ex K_, GiNaC::ex ws_, GiNaC_symbol_set lm_,
                             GiNaC_symbol_set em_, subs_list rm_, symbol_factory& sf_)
  : tm(std::move(tm_)),
    K(std::move(K_)),
    WickProduct(std::move(ws_)),
    loop_momenta(std::move(lm_)),
    external_momenta(std::move(em_)),
    Rayleigh_momenta(std::move(rm_)),
    sf(sf_)
  {
    // no need to apply any transformations if this is a tree-level term
    if(loop_momenta.empty()) return;

    // STEP 1 - pair nontrivial arguments appearing in the Wick product with Rayleigh momenta.
    // In general, we want integrals over the power spectra to have a simple argument that
    // corresponds to the magnitude of one of the integration variables
    this->match_Wick_to_Rayleigh();

    // STEP 2 - canonicalize all labels
    // This relabels all loop momenta and Rayleigh momenta to standard form.
    this->canonicalize_loop_labels();
    this->canonicalize_Rayleigh_labels();
  }


void loop_integral::write(std::ostream& out) const
  {
    std::cout << "  time function = " << this->tm << '\n';
    std::cout << "  momentum kernel = " << this->K << '\n';
    std::cout << "  Wick product = " << this->WickProduct << '\n';

    if(!this->loop_momenta.empty())
      {
        std::cout << "  loop momenta =";
        for(const auto& sym : this->loop_momenta)
          {
            std::cout << " " << sym;
          }
        std::cout << '\n';
      }

    if(!this->Rayleigh_momenta.empty())
      {
        std::cout << "  Rayleigh momenta =";
        for(const auto& rule : this->Rayleigh_momenta)
          {
            std::cout << " " << rule.first << " -> " << rule.second << ";";
          }
        std::cout << '\n';
      }
  }


void loop_integral::canonicalize_loop_labels()
  {
    GiNaC_symbol_set new_momenta;
    GiNaC::exmap relabel;

    // step through loop momenta, relabelling to canonicalized variables
    unsigned int count = 0;
    for(const auto& l : this->loop_momenta)
      {
        const auto L = this->sf.make_canonical_loop_momentum(count++);
        relabel[l] = L;
        new_momenta.insert(L);
      }

    // propagate this relabelling to the kernel and the Wick product
    this->K = this->K.subs(relabel);
    this->WickProduct = this->WickProduct.subs(relabel);

    // propagate to any Rayleigh replacement rules in use
    for(auto& rule : this->Rayleigh_momenta)
      {
        rule.second = rule.second.subs(relabel);
      }

    this->loop_momenta = new_momenta;
  }


void loop_integral::canonicalize_Rayleigh_labels()
  {
    GiNaC::exmap new_Rayleigh;
    GiNaC::exmap relabel;

    // step through Rayleigh momenta, relabelling to canonicalized forms
    unsigned int count = 0;
    for(const auto& rule : this->Rayleigh_momenta)
      {
        const auto S = this->sf.make_canonical_Rayleigh_momentum(count++);
        relabel[rule.first] = S;
        new_Rayleigh[S] = rule.second;    // RHS of rules never depend on other Rayleigh momenta, so this is safe
      }

    // propagate this relabelling to the kernel and Wick product
    this->K = this->K.subs(relabel);
    this->WickProduct = this->WickProduct.subs(relabel);

    this->Rayleigh_momenta = new_Rayleigh;
  }


void loop_integral::match_Wick_to_Rayleigh()
  {
    GiNaC::ex new_Wick{1};

    if(!GiNaC::is_a<GiNaC::mul>(this->WickProduct))
      throw exception(ERROR_BADLY_FORMED_WICK_PRODUCT, exception_code::loop_transformation_error);

    for(size_t i = 0; i < this->WickProduct.nops(); ++i)
      {
        const auto& factor = this->WickProduct.op(i);
        const auto& f = GiNaC::ex_to<GiNaC::function>(factor);

        if(f.get_name() == "Pk")
          {
            const auto& f1 = GiNaC::ex_to<GiNaC::symbol>(f.op(0));
            const auto& f2 = GiNaC::ex_to<GiNaC::symbol>(f.op(1));
            const auto& arg = f.op(2);

            // if argument is a simple symbol (a loop momentum or a Rayleigh momentum, or an external momentum) then we have nothing to do
            if(GiNaC::is_a<GiNaC::symbol>(arg))
              {
                const auto& sym = GiNaC::ex_to<GiNaC::symbol>(arg);
                if(this->loop_momenta.find(sym) != this->loop_momenta.end()) { new_Wick *= cfs::Pk(f1, f2, sym); continue; }
                if(this->Rayleigh_momenta.find(sym) != this->Rayleigh_momenta.end()) { new_Wick *= cfs::Pk(f1, f2, sym); continue; }
                if(this->external_momenta.find(sym) != this->external_momenta.end()) { new_Wick *= cfs::Pk(f1, f2, sym); continue; }

                std::cerr << factor << '\n';

                throw exception(ERROR_UNKNOWN_WICK_PRODUCT_LABEL, exception_code::loop_transformation_error);
              }

            // otherwise, try to match this argument to a Rayleigh momentum (or its negative)
            auto t = this->Rayleigh_momenta.begin();
            for(; t != this->Rayleigh_momenta.end(); ++t)
              {
                // if a direct match, just relabel the argument; don't need to worry about sign inversion
                // since the power spectrum just depends on |arg|
                if(static_cast<bool>(t->second == arg)
                   || static_cast<bool>(t->second == -arg))
                  {
                    new_Wick *= cfs::Pk(f1, f2, t->first);
                    break;
                  }
              }

            if(t == this->Rayleigh_momenta.end())
              {
                if(!this->Rayleigh_momenta.empty())
                  {
                    std::cerr << factor << '\n';
                    throw exception(ERROR_CANT_MATCH_WICK_TO_RAYLEIGH, exception_code::loop_transformation_error);
                  }

                // matching Rayleigh momentum didn't exist, so insert one
                auto S = this->sf.make_canonical_Rayleigh_momentum(0);
                new_Wick *= cfs::Pk(f1, f2, S);
                this->Rayleigh_momenta[S] = arg;
              }

            continue;
          }

        throw exception(ERROR_BADLY_FORMED_WICK_PRODUCT, exception_code::loop_transformation_error);
      }

    // replace old Wick product with matched version
    this->WickProduct = new_Wick;
  }
