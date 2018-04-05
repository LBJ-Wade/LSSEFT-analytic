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


#include <string>
#include <sstream>

#include "GiNaC_print.h"

#include "localizations/messages.h"
#include "shared/exceptions.h"


// forward-declare print functions
static std::string print_func(const GiNaC::ex& expr, const std::string& name);
static std::string print_ginac(const GiNaC::ex& expr);
static std::string print_operands(const GiNaC::ex& expr, const std::string& op);
static std::string print_constant(const GiNaC::ex& expr);


static const std::map< std::string, std::string > func_convert
  {
    {"abs", "Abs"},
    {"sqrt", "Sqrt"},
    {"sin", "Sin"},
    {"cos", "Cos"},
    {"tan", "Tan"},
    {"asin", "ArcSin"},
    {"acos", "ArcCos"},
    {"atan", "ArcTan"},
    {"sinh", "Sinh"},
    {"cosh", "Cosh"},
    {"tanh", "Tanh"},
    {"asinh", "ArcSinh"},
    {"acosh", "ArcCosh"},
    {"atanh", "ArcTanh"},
    {"exp", "Exp"},
    {"log", "Log"},
    {"power", "Power"},
    {"FabJ", "FabJ"}
  };


static std::string print_func(const GiNaC::ex& expr, const std::string& name)
  {
    auto t = func_convert.find(name);
    if(t == func_convert.end())
      {
        std::ostringstream msg;
        msg << ERROR_UNKNOWN_GINAC_FUNCTION << " '" << name << "'";
        throw exception(msg.str(), exception_code::backend_error);
      }

    std::string rval{t->second};
    rval.append("[");
    rval.append(print_operands(expr, ","));
    rval.append("]");

    return rval;
  }


static std::string print_ginac(const GiNaC::ex& expr)
  {
    std::ostringstream out;
    out << expr;
    return out.str();
  }


static std::string print_operands(const GiNaC::ex& expr, const std::string& op)
  {
    std::string rval;

    unsigned int c = 0;
    for(const auto& arg : expr)
      {
        if(c > 0) rval.append(op);

        // determine whether we should bracket this operand
        // this should happen any time the operand has lower precedence than the current operand, but at the moment
        // we only deal with + and * so we can simply check for +
        // we can kill the brackets if the operands are a comma-separated list of arguments, though
        bool bracket = op != "," && GiNaC::is_a<GiNaC::add>(arg);

        if(bracket) rval.append("(");
        rval.append(format_print(arg));
        if(bracket) rval.append(")");

        ++c;
      }

    return rval;
  }


static std::string print_constant(const GiNaC::ex& expr)
  {
    const auto& c = GiNaC::ex_to<GiNaC::constant>(expr);

    std::ostringstream buf;
    buf << c;

    if(buf.str() == "Pi") return std::string{"Pi"};
    return buf.str();
  }


std::string format_print(const GiNaC::ex& expr)
  {
    std::string name;

    if(GiNaC::is_a<GiNaC::function>(expr)) name = GiNaC::ex_to<GiNaC::function>(expr).get_name();
    else name = GiNaC::ex_to<GiNaC::basic>(expr).class_name();

    if     (name == "numeric")   return print_ginac(expr);
    else if(name == "symbol")    return print_ginac(expr);
    else if(name == "add")       return print_operands(expr, "+");
    else if(name == "mul")       return print_operands(expr, "*");
    else if(name == "constant")  return print_constant(expr);
    else if(name == "tensdelta") return print_ginac(expr);
    else if(name == "idx")       return print_ginac(expr);
    else if(name == "varidx")    return print_ginac(expr);
    else if(name == "indexed")   return print_ginac(expr);
    else if(name == "D")         return std::string{"Dz"};
    else if(name == "DA")        return std::string{"DAz"};
    else if(name == "DB")        return std::string{"DBz"};
    else if(name == "DD")        return std::string{"DDz"};
    else if(name == "DE")        return std::string{"DEz"};
    else if(name == "DF")        return std::string{"DFz"};
    else if(name == "DG")        return std::string{"DGz"};
    else if(name == "DJ")        return std::string{"DJx"};
    else if(name == "f")         return std::string("fz");
    else if(name == "fA")        return std::string("fAz");
    else if(name == "fB")        return std::string("fBz");
    else if(name == "fD")        return std::string("fDz");
    else if(name == "fE")        return std::string("fEz");
    else if(name == "fF")        return std::string("fFz");
    else if(name == "fG")        return std::string("fGz");
    else if(name == "fJ")        return std::string("fJz");

    // not a standard operation, so assume it must be a special function
    // look up its C++ form in func_map, and then format its arguments,
    // taking care to keep track of use counts
    return print_func(expr, name);
  }
