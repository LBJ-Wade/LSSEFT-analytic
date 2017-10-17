//
// Created by David Seery on 17/10/2017.
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

#include <iomanip>
#include <sstream>

#include "formatter.h"
#include "localizations/messages.h"


std::string format_time(boost::timer::nanosecond_type time, unsigned int precision)
  {
    std::ostringstream out;

    constexpr boost::timer::nanosecond_type mu_sec = 1000;
    constexpr boost::timer::nanosecond_type m_sec  = 1000*mu_sec;
    constexpr boost::timer::nanosecond_type sec    = 1000*m_sec;
    constexpr boost::timer::nanosecond_type minute = 60*sec;
    constexpr boost::timer::nanosecond_type hour   = 60*minute;

    if(time > hour)
      {
        out << time/hour << FORMAT_HOUR_LABEL << " ";
        time = time % hour;
      }
    if(time > minute)
      {
        out << time/minute << FORMAT_MINUTE_LABEL << " ";
        time = time % minute;
      }
    out << std::setprecision(precision) << static_cast<double>(time) / sec << FORMAT_SECOND_LABEL;

    return(out.str());
  }
