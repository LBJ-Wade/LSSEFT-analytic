//
// Created by David Seery on 23/08/2017.
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

#ifndef LSSEFT_ANALYTIC_ERROR_H
#define LSSEFT_ANALYTIC_ERROR_H


#include <string>


class error_handler
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor is default
    error_handler() = default;
    
    //! destructor is default
    ~error_handler() = default;
    
    
    // INTERFACE
  
  public:
    
    //! report an error
    void error(std::string msg);
    
    //! report a warning
    void warn(std::string msg);
    
    //! report an information message
    void info(std::string msg);
    
    //! make an announcement
    void announce(std::string msg);
    
  };


#endif //LSSEFT_ANALYTIC_ERROR_H
