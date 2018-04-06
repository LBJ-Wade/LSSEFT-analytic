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

#ifndef LSSEFT_ANALYTIC_RAYLEIGH_MOMENTA_H
#define LSSEFT_ANALYTIC_RAYLEIGH_MOMENTA_H


#include "services/service_locator.h"
#include "utilities/GiNaC_utils.h"


namespace Rayleigh
  {
    
    //! merge substitution rules from 'source' into an existing map 'dest', applying substitutions in
    //! 'subs_rules' to the RHS.
    //! returns a GiNaC::exmap which performs any relabelling needed due to symbol collisions
    //! with symbols already in map, or in the set 'reserved_symbols'
    GiNaC::exmap
    merge_Rayleigh_rules(GiNaC::exmap& dest, const GiNaC::exmap& source, const GiNaC_symbol_set& reserved,
                         const GiNaC::exmap& subs_rules, service_locator& loc);
    
    //! convenience wrapper to drive merge_Rayleigh_rules(), updating reserved symbols as we go
    GiNaC::exmap merge_Rayleigh_lists(const GiNaC::exmap& source, GiNaC::exmap& dest, GiNaC_symbol_set& reserved,
                                      const GiNaC::exmap& subs_rules, service_locator& loc);
    
    //! build a replacement rule list for trivial Rayleigh labels
    GiNaC::exmap remove_Rayleigh_trivial(GiNaC::exmap& rm);
    
    //! prune a list of Rayleigh labels
    void prune_Rayleigh_list(GiNaC::exmap& list, const GiNaC::ex& K);
    
  }   // namespace Rayleigh


#endif //LSSEFT_ANALYTIC_RAYLEIGH_MOMENTA_H
