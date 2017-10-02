//
// Created by David Seery on 01/10/2017.
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


#include <iostream>

#include "argument_cache.h"
#include "switches.h"

#include "shared/common.h"

#include "boost/program_options.hpp"


argument_cache::argument_cache(int argc, char**& argv)
  {
    // set up BOOST::program_options descriptors for command-line arguments
    boost::program_options::options_description generic{"Generic"};
    generic.add_options()
      (SWITCH_HELP, HELP_HELP)
      (SWITCH_VERSION, HELP_VERSION)
      ;

    boost::program_options::options_description expressions{"Expression handling"};
    expressions.add_options()
      (SWITCH_AUTO_SYMMETRIZE, HELP_AUTO_SYMMETRIZE)
      (SWITCH_22_SYMMETRIZE, HELP_22_SYMMETRIZE)
      ;

    boost::program_options::options_description backend{"Backend control"};
    backend.add_options()
      (SWITCH_COUNTERTERMS, HELP_COUNTERTERMS)
      (SWITCH_OUTPUT, boost::program_options::value<std::string>(), HELP_OUTPUT)
      ;

    boost::program_options::options_description backend_hidden{"Hidden backed control options"};
    backend_hidden.add_options()
      (SWITCH_NO_COUNTERTERMS, "")
      ;

    boost::program_options::options_description expressions_hidden{"Hidden expression options"};
    expressions_hidden.add_options()
      (SWITCH_NO_AUTO_SYMMETRIZE, "")
      (SWITCH_NO_22_SYMMETRIZE, "")
      ;

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(expressions).add(expressions_hidden).add(backend);

    boost::program_options::options_description output_options;
    output_options.add(generic).add(expressions).add(backend);

    boost::program_options::variables_map option_map;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmdline_options), option_map);
    boost::program_options::notify(option_map);

    bool emitted_version = false;

    if(option_map.count(SWITCH_VERSION))
      {
        std::cout << PROGRAM_NAME << " " << PROGRAM_VERSION << " " << PROGRAM_COPYRIGHT << '\n';
        emitted_version = true;
      }

    if(option_map.count(SWITCH_HELP))
      {
        if(!emitted_version) std::cout << PROGRAM_NAME << " " << PROGRAM_VERSION << " " << PROGRAM_COPYRIGHT << '\n';
        std::cout << output_options << '\n';
        exit(EXIT_SUCCESS);
      }

    if(option_map.count(SWITCH_AUTO_SYMMETRIZE))    this->auto_symmetrize = true;
    if(option_map.count(SWITCH_NO_AUTO_SYMMETRIZE)) this->auto_symmetrize = false;
    if(option_map.count(SWITCH_22_SYMMETRIZE))      this->symmetrize_22 = true;
    if(option_map.count(SWITCH_NO_22_SYMMETRIZE))   this->symmetrize_22 = false;

    if(option_map.count(SWITCH_COUNTERTERMS))       this->counterterms = true;
    if(option_map.count(SWITCH_NO_COUNTERTERMS))    this->counterterms = false;

    if(option_map.count(SWITCH_OUTPUT_LONG))
      {
        boost::filesystem::path outpath = option_map[SWITCH_OUTPUT_LONG].as<std::string>();
        if(!outpath.is_absolute()) outpath = boost::filesystem::absolute(outpath);

        this->output_root = std::move(outpath);
      }
  }


bool argument_cache::get_auto_symmetrize() const
  {
    return this->auto_symmetrize;
  }


bool argument_cache::get_symmetrize_22() const
  {
    return this->symmetrize_22;
  }


const boost::filesystem::path& argument_cache::get_output_root() const
  {
    return this->output_root;
  }


bool argument_cache::get_counterterms() const
  {
    return this->counterterms;
  }


