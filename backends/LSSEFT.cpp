//
// Created by David Seery on 02/10/2017.
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
#include <fstream>
#include <sstream>

#include "LSSEFT.h"

#include "utilities/GiNaC_utils.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


namespace LSSEFT_impl
  {

    LSSEFT_kernel::LSSEFT_kernel(GiNaC::ex ig_, GiNaC::ex ms_, GiNaC::ex wp_, GiNaC_symbol_set vs_, GiNaC_symbol_set em_)
      : integrand(std::move(ig_)),
        measure(std::move(ms_)),
        WickProduct(std::move(wp_)),
        variables(std::move(vs_)),
        external_momenta(std::move(em_))
      {
      }


    bool LSSEFT_kernel::is_equal(const LSSEFT_kernel& obj) const
      {
        // coalesce measure and integrand, then test for equality
        auto combined_a = (this->measure * this->integrand).expand();
        auto combined_b = (obj.measure * obj.integrand).expand();

        if(!static_cast<bool>(combined_a == combined_b)) return false;

        // test for equality of Wick product
        if(!static_cast<bool>(this->WickProduct == obj.WickProduct)) return false;


        // test for equality of integration variables
        auto a_iv = order_symbol_set(this->variables);
        auto b_iv = order_symbol_set(obj.variables);

        if(!std::equal(a_iv.cbegin(), a_iv.cend(), b_iv.cbegin(), b_iv.cend(),
                       [](const GiNaC::symbol& asym, const GiNaC::symbol& bsym) -> bool
                         { return asym.get_name() == bsym.get_name(); })) return false;

        // test for equality of external momenta
        auto a_em = order_symbol_set(this->external_momenta);
        auto b_em = order_symbol_set(obj.external_momenta);

        if(!std::equal(a_em.cbegin(), a_em.cend(), b_em.cbegin(), b_em.cend(),
                       [](const GiNaC::symbol& asym, const GiNaC::symbol& bsym) -> bool
                         { return asym.get_name() == bsym.get_name(); })) return false;

        return true;
      }


    size_t LSSEFT_kernel::hash() const
      {
        size_t h = 0;

        // coalesce measure and integrand, then expand to get in a canonical form
        auto combined = (this->measure * this->integrand).expand();

        // print to string and hash
        std::ostringstream expr_string;
        expr_string << combined;

        hash_impl::hash_combine(h, expr_string.str());

        // order integration variables lexically, convert to a string, and hash
        auto ordered_iv = order_symbol_set(this->variables);

        std::string iv_string;
        std::for_each(ordered_iv.begin(), ordered_iv.end(),
                      [&](const GiNaC::symbol& e) -> void
                        { iv_string += e.get_name(); });

        hash_impl::hash_combine(h, iv_string);

        // order external momenta lexically, convert to a string, and hash
        auto ordered_em = order_symbol_set(this->external_momenta);

        std::string em_string;
        std::for_each(ordered_em.begin(), ordered_em.end(),
                      [&](const GiNaC::symbol& e) -> void
                        { em_string += e.get_name(); });

        hash_impl::hash_combine(h, em_string);

        return h;
      }


    static const std::map< std::string, std::string > func_convert
      {
        {"abs", "std::abs"},
        {"sqrt", "std::sqrt"},
        {"sin", "std::sin"},
        {"cos", "std::cos"},
        {"tan", "std::tan"},
        {"asin", "std::asin"},
        {"acos", "std::acos"},
        {"atan", "std::atan"},
        {"atan2", "std::atan2"},
        {"sinh", "std::sinh"},
        {"cosh", "std::cosh"},
        {"tanh", "std::tanh"},
        {"asinh", "std::asinh"},
        {"acosh", "std::acosh"},
        {"atanh", "std::atanh"},
        {"exp", "std::exp"},
        {"log", "std::log"},
        {"pow", "std::pow"},
        {"tgamma", "std::tgamma"},
        {"lgamma", "std::lgamma"}
      };


    // forward-declare print function
    static std::string format_print(const GiNaC::ex& expr);


    static std::string print_ginac(const GiNaC::ex& expr)
      {
        std::ostringstream out;
        out << GiNaC::csrc << expr;
        return out.str();
      }


    static std::string print_operands(const GiNaC::ex& expr, const std::string& op)
      {
        std::string rval;

        unsigned int c = 0;
        for(const auto& arg : expr)
          {
            if(c > 0) rval.append(op);
            rval.append(format_print(arg));

            ++c;
          }

        return rval;
      }


    // treat powers specially, because it's better to unroll them and
    // avoid a function call if we have an exact integral power
    static std::string print_power(const GiNaC::ex& expr)
      {
        std::ostringstream out;
        size_t n = expr.nops();

        if(n != 2) throw exception(ERROR_BACKEND_POW_ARGUMENTS, exception_code::backend_error);

        // perform use-counting on exponent, which is necessary since its values may have been
        // used elsewhere in the CSE tree, even if it is an integer that we will unroll and not use explicitly
        std::string exponent = format_print(expr.op(1));

        const GiNaC::ex& exp_generic = expr.op(1);

        if(GiNaC::is_a<GiNaC::numeric>(exp_generic))
          {
            const auto& exp_numeric = GiNaC::ex_to<GiNaC::numeric>(exp_generic);

            std::string base = format_print(expr.op(0));

            if(GiNaC::is_integer(exp_numeric))
              {

                if(GiNaC::is_nonneg_integer(exp_numeric))
                  {
                    if(exp_numeric.to_int() == 0) out << "1.0";
                    else if(exp_numeric.to_int() == 1) out << base;
                    else if(exp_numeric.to_int() == 2) out << base << "*" << base;
                    else if(exp_numeric.to_int() == 3) out << base << "*" << base << "*" << base;
                    else if(exp_numeric.to_int() == 4) out << base << "*" << base << "*" << base << "*" << base;
                    else out << "std::pow(" << base << "," << exp_numeric.to_int() << ")";
                  }
                else  // negative integer
                  {
                    out << "1.0/";
                    if(exp_numeric.to_int() == -0) out << "1.0";
                    else if(exp_numeric.to_int() == -1) out << base;
                    else if(exp_numeric.to_int() == -2) out << "(" << base << "*" << base << ")";
                    else if(exp_numeric.to_int() == -3) out << "(" << base << "*" << base << "*" << base << ")";
                    else if(exp_numeric.to_int() == -4)
                      out << "(" << base << "*" << base << "*" << base << "*" << base << ")";
                    else out << "std::pow(" << base << "," << -exp_numeric.to_int() << ")";
                  }
              }
            else  // not an integer
              {
                out << "std::pow(" << base << "," << exponent << ")";
              }
          }

        return(out.str());
      }


    static std::string print_Pk(const GiNaC::ex& expr)
      {
        std::ostringstream out;
        size_t n = expr.nops();

        if(n != 3) throw exception(ERROR_BACKEND_PK_ARGUMENTS, exception_code::backend_error);

        out << "data_->Pk(data_->IR_cutoff, data_->UV_cutoff, " << format_print(expr.op(2)) << ")";

        return out.str();
      }


    static std::string format_print(const GiNaC::ex& expr)
      {
        std::string name;

        if(GiNaC::is_a<GiNaC::function>(expr)) name = GiNaC::ex_to<GiNaC::function>(expr).get_name();
        else name = GiNaC::ex_to<GiNaC::basic>(expr).class_name();

        if     (name == "numeric")   return print_ginac(expr);
        else if(name == "symbol")    return print_ginac(expr);
        else if(name == "add")       return print_operands(expr, "+");
        else if(name == "mul")       return print_operands(expr, "*");
        else if(name == "power")     return print_power(expr);
        else if(name == "constant")  return print_ginac(expr);
        else if(name == "tensdelta") return print_ginac(expr);
        else if(name == "idx")       return print_ginac(expr);
        else if(name == "varidx")    return print_ginac(expr);
        else if(name == "indexed")   return print_ginac(expr);

        // not a standard operation, so assume it must be a special function
        // look up its C++ form in func_map, and then format its arguments,
        // taking care to keep track of use counts

        if(name == "Pk") return print_Pk(expr);

        auto t = func_convert.find(name);
        if(t == func_convert.end())
          {
            std::ostringstream msg;
            msg << ERROR_UNKNOWN_GINAC_FUNCTION << " '" << name << "'";
            throw exception(msg.str(), exception_code::backend_error);
          }

        std::string rval{t->second};
        rval.append("(");
        rval.append(print_operands(expr, ","));
        rval.append(")");

        return rval;
      }


    std::string LSSEFT_kernel::print_integrand(const GiNaC::exmap& subs_map) const
      {
        return format_print(this->integrand.subs(subs_map));
      }


    std::string LSSEFT_kernel::print_measure(const GiNaC::exmap& subs_map) const
      {
        return format_print(this->measure.subs(subs_map));
      }


    std::string LSSEFT_kernel::print_WickProduct(const GiNaC::exmap& subs_map) const
      {
        return format_print(this->WickProduct.subs(subs_map));
      }

  }   // namespace LSSEFT_impl


LSSEFT::LSSEFT(boost::filesystem::path rt_, service_locator& lc_)
  : loc(lc_),
    root(std::move(rt_))
  {
  }


LSSEFT& LSSEFT::add(const Pk_rsd& P, std::string name)
  {
    // check whether a power spectrum with this name has already been registerd
    auto it = this->Pk_db.find(name);

    if(it != this->Pk_db.end())  // trouble -- one already exists
      {
        std::ostringstream msg;
        msg << ERROR_BACKEND_PK_RSD_ALREADY_REGISTERED_A
            << " '" << name << "' "
            << ERROR_BACKEND_PK_RSD_ALREADY_REGISTERED_B;
        throw exception(msg.str(), exception_code::backend_error);
      }

    // no power spectrum with this name was found, so insert a new database record with a
    // reference to P
    auto res = this->Pk_db.insert(std::make_pair(std::move(name), std::cref(P)));
    if(!res.second) throw exception(ERROR_BACKEND_PK_INSERT_FAILED, exception_code::backend_error);

    // process P's kernels
    this->process_kernels(P.get_13());
    this->process_kernels(P.get_22());

    return *this;
  }


void LSSEFT::process_kernels(const Pk_rsd_group& group)
  {
    auto visitor = [&](const one_loop_element& elt) -> void
      {
        using LSSEFT_impl::LSSEFT_kernel;

        LSSEFT_kernel ker{elt.get_integrand(), elt.get_measure(), elt.get_Wick_product(),
                          elt.get_integration_variables(), elt.get_external_momenta()};

        auto it = this->kernel_db.find(ker);

        if(it != this->kernel_db.end()) return;

        // generate a new kernel name
        auto res = this->kernel_db.insert(std::make_pair(std::move(ker), this->make_unique_kernel_name()));
        if(!res.second) throw exception(ERROR_BACKEND_KERNEL_INSERT_FAILED, exception_code::backend_error);
      };

    // visit the kernel elements associated with each power of mu
    // and insert them into our kernel database
    group.visit({0,2,4,6,8}, visitor);
  }


void LSSEFT::write() const
  {
    // generate create statements
    this->write_create();

    // generate kernel integrands
    this->write_kernels();
  }


boost::filesystem::path LSSEFT::make_output_path(std::string append) const
  {
    return this->root.parent_path() / boost::filesystem::path{this->root.stem().string() + "_" + append + ".cpp"};
  }


void LSSEFT::write_create() const
  {
    auto output = this->make_output_path("create");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->kernel_db)
      {
        const std::string& name = record.second;

        // write create statements for all kernels that we require
        outf << "create_impl::oneloop_momentum_integral_table(db, \"" << name << "\", policy);" << '\n';
      }
    outf << '\n';

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        // first, write create statements for the power spectrum mu coefficients
        outf << "create_impl::oneloop_rsd_Pk_table(db, \"" << name << "_mu0\", policy);" << '\n';
        outf << "create_impl::oneloop_rsd_Pk_table(db, \"" << name << "_mu2\", policy);" << '\n';
        outf << "create_impl::oneloop_rsd_Pk_table(db, \"" << name << "_mu4\", policy);" << '\n';
        outf << "create_impl::oneloop_rsd_Pk_table(db, \"" << name << "_mu6\", policy);" << '\n';
        outf << "create_impl::oneloop_rsd_Pk_table(db, \"" << name << "_mu8\", policy);" << '\n';
        outf << '\n';

        // second, the projection to P0, P2, P4
        outf << "create_impl::multipole_Pk_table(db, \"" << name << "_P0\", policy);" << '\n';
        outf << "create_impl::multipole_Pk_table(db, \"" << name << "_P2\", policy);" << '\n';
        outf << "create_impl::multipole_Pk_table(db, \"" << name << "_P4\", policy);" << '\n';
        outf << '\n';
      }

    outf.close();
  }


void LSSEFT::write_kernels() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("kernels");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    auto& sf = this->loc.get_symbol_factory();
    auto q0 = sf.make_canonical_loop_momentum(0);
    auto x = sf.make_symbol("x");

    auto q_ = sf.make_symbol("q_");
    auto z_ = sf.make_symbol("z_");
    auto k_ = sf.make_symbol("k_");

    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        const auto& external_momenta = kernel.get_external_momenta();
        const auto& k = *external_momenta.begin();

        GiNaC::exmap subs_map = { {q0, q_}, {x, z_}, {k, k_} };

        // write create statements for all kernels that we require
        outf << "static int " << name << "_integrand(const int* ndim_, const cubareal x_[], const int* ncomp_, cubareal f_[], void* userdata_)" << '\n';
        outf << " {" << '\n';
        outf << "   oneloop_momentum_impl::integrand_data* data_ = static_cast<integrand_data*>(userdata_);" << '\n';
        outf << '\n';
        outf << "   Mpc_units::energy k_ = data_->k;" << '\n';
        outf << "   Mpc_units::energy q_ = data_->IR_cutoff + x_[0] * data_->q_range;" << '\n';
        outf << "   double z_ = 2.0*x_[1] - 1.0;" << '\n';
        outf << '\n';
        outf << "   auto value_ = " << kernel.print_integrand(subs_map) << ";" << '\n';
        outf << "   value_ *= " << kernel.print_measure(subs_map) << ";" << '\n';
        outf << "   value_ *= " << kernel.print_WickProduct(subs_map) << ";" << '\n';
        outf << "   f_[0] = data_->jacobian_d3q * value_;" << '\n';
        outf << '\n';
        outf << "   return 0;" << '\n';
        outf << " }" << '\n';
      }

    outf.close();
  }


std::string LSSEFT::make_unique_kernel_name()
  {
    return this->kernel_root + std::to_string(this->kernel_count++);
  }
