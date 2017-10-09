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

    LSSEFT_kernel::LSSEFT_kernel(GiNaC::ex ig_, GiNaC::ex ms_, GiNaC::ex wp_, GiNaC_symbol_set vs_, GiNaC_symbol_set em_,
                                     mass_dimension dm_)
      : integrand(std::move(ig_)),
        measure(std::move(ms_)),
        WickProduct(std::move(wp_)),
        variables(std::move(vs_)),
        external_momenta(std::move(em_)),
        dim(dm_)
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


    // forward-declare print functions
    static std::string format_print(const GiNaC::ex& expr);
    static std::string print_operands(const GiNaC::ex& expr, const std::string& op);


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
        rval.append("(");
        rval.append(print_operands(expr, ","));
        rval.append(")");

        return rval;
      }


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


    static std::string unwrap_power(const GiNaC::ex& base, unsigned int factors)
      {
        bool simple = GiNaC::is_a<GiNaC::symbol>(base) || GiNaC::is_a<GiNaC::constant>(base);
        std::ostringstream str;

        std::string base_str = format_print(base);

        unsigned int count = 0;
        while(factors > 0)
          {
            if(count > 0) str << "*";
            if(simple) str << base_str;
            else       str << "(" << base_str << ")";

            --factors;
            ++count;
          }

        return str.str();
      }


    // treat powers specially, because it's better to unroll them and
    // avoid a function call if we have an exact integral power
    static std::string print_power(const GiNaC::ex& expr)
      {
        std::ostringstream out;
        size_t n = expr.nops();

        if(n != 2) throw exception(ERROR_BACKEND_POW_ARGUMENTS, exception_code::backend_error);

        const GiNaC::ex& base_expr = expr.op(0);
        const GiNaC::ex& exp_expr = expr.op(1);

        if(GiNaC::is_a<GiNaC::numeric>(exp_expr))
          {
            const auto& exp_numeric = GiNaC::ex_to<GiNaC::numeric>(exp_expr);

            if(GiNaC::is_integer(exp_numeric))
              {
                if(GiNaC::is_nonneg_integer(exp_numeric))
                  {
                    if(exp_numeric.to_int() == 0) out << "1.0";
                    else if(exp_numeric.to_int() == 1) out << unwrap_power(base_expr, 1);
                    else if(exp_numeric.to_int() == 2) out << "(" << unwrap_power(base_expr, 2) << ")";
                    else if(exp_numeric.to_int() == 3) out << "(" << unwrap_power(base_expr, 3) << ")";
                    else if(exp_numeric.to_int() == 4) out << "(" << unwrap_power(base_expr, 4) << ")";
                    else out << "std::pow(" << format_print(base_expr) << "," << exp_numeric.to_int() << ")";
                  }
                else  // negative integer
                  {
                    if(exp_numeric.to_int() == -0) out << "1.0";
                    else if(exp_numeric.to_int() == -1) out << "1.0/" << unwrap_power(base_expr, 1);
                    else if(exp_numeric.to_int() == -2) out << "1.0/" << "(" << unwrap_power(base_expr, 2) << ")";
                    else if(exp_numeric.to_int() == -3) out << "1.0/" << "(" << unwrap_power(base_expr, 3) << ")";
                    else if(exp_numeric.to_int() == -4) out << "1.0/" << "(" << unwrap_power(base_expr, 4) << ")";
                    else out << "std::pow(" << format_print(base_expr) << "," << exp_numeric.to_int() << ")";
                  }
              }
            else  // not an integer
              {
                out << "std::pow(" << format_print(base_expr) << "," << format_print(exp_expr) << ")";
              }
          }

        return(out.str());
      }


    static std::string print_Pk(const GiNaC::ex& expr)
      {
        std::ostringstream out;
        size_t n = expr.nops();

        if(n != 3) throw exception(ERROR_BACKEND_PK_ARGUMENTS, exception_code::backend_error);

        out << "(data_->Pk((" << format_print(expr.op(2)) << ")/Mpc_units::Mpc) / Mpc_units::Mpc3)";

        return out.str();
      }


    static std::string print_constant(const GiNaC::ex& expr)
      {
        const auto& c = GiNaC::ex_to<GiNaC::constant>(expr);

        std::ostringstream buf;
        buf << c;

        if(buf.str() == "Pi") return std::string{"M_PI"};
        return buf.str();
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
        else if(name == "constant")  return print_constant(expr);
        else if(name == "tensdelta") return print_ginac(expr);
        else if(name == "idx")       return print_ginac(expr);
        else if(name == "varidx")    return print_ginac(expr);
        else if(name == "indexed")   return print_ginac(expr);
        else if(name == "Pk")        return print_Pk(expr);
        else if(name == "D")         return std::string{"val.D_lin"};
        else if(name == "DA")        return std::string{"val.A"};
        else if(name == "DB")        return std::string{"val.B"};
        else if(name == "DD")        return std::string{"val.D"};
        else if(name == "DE")        return std::string{"val.E"};
        else if(name == "DF")        return std::string{"val.F"};
        else if(name == "DG")        return std::string{"val.G"};
        else if(name == "DJ")        return std::string{"val.J"};
        else if(name == "f")         return std::string("val.f_lin");
        else if(name == "fA")        return std::string("val.fA");
        else if(name == "fB")        return std::string("val.fB");
        else if(name == "fD")        return std::string("val.fD");
        else if(name == "fE")        return std::string("val.fE");
        else if(name == "fF")        return std::string("val.fF");
        else if(name == "fG")        return std::string("val.fG");
        else if(name == "fJ")        return std::string("val.fJ");

        // not a standard operation, so assume it must be a special function
        // look up its C++ form in func_map, and then format its arguments,
        // taking care to keep track of use counts
        return print_func(expr, name);
      }


    std::string LSSEFT_kernel::print_integrand(const GiNaC::exmap& subs_map) const
      {
        return format_print(this->integrand.subs(subs_map));
      }


    std::string LSSEFT_kernel::print_measure(const GiNaC::exmap& subs_map) const
      {
        // LSSEFT uses a constant normalization of 8pi^2 for the integrals
        return format_print(this->measure.subs(subs_map));
      }


    GiNaC::ex filter(const GiNaC::ex& factor, const GiNaC_symbol_set& external_momenta)
      {
        const auto& func = GiNaC::ex_to<GiNaC::function>(factor);

        if(func.get_name() == "Pk" && func.nops() == 3)
          {
            const GiNaC::ex& arg = func.op(2);
            if(GiNaC::is_a<GiNaC::symbol>(arg))
              {
                const auto& sym = GiNaC::ex_to<GiNaC::symbol>(arg);
                auto it = external_momenta.find(sym);

                // if symbol included in external momenta, return unity
                if(it != external_momenta.end()) return GiNaC::ex{1};
              }
          }

        return factor;
      }


    std::string
    LSSEFT_kernel::print_WickProduct(const GiNaC::exmap& subs_map, const GiNaC_symbol_set& external_momenta) const
      {
        GiNaC::ex filtered{1};

        // filter out factors in the Wick product that don't participate in the integration
        if(GiNaC::is_a<GiNaC::function>(this->WickProduct))
          {
            filtered *= filter(this->WickProduct, external_momenta);
          }
        else if(GiNaC::is_a<GiNaC::mul>(this->WickProduct))
          {
            size_t ops = this->WickProduct.nops();
            for(size_t i = 0; i < ops; ++i)
              {
                const GiNaC::ex& factor = this->WickProduct.op(i);

                if(GiNaC::is_a<GiNaC::function>(factor))
                  {
                    filtered *= filter(factor, external_momenta);
                  }
                else
                  {
                    filtered *= factor;
                  }
              }
          }

        return format_print(filtered.subs(subs_map));
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
    using LSSEFT_impl::mass_dimension;
    this->process_kernels(P.get_13(), mass_dimension::zero);
    this->process_kernels(P.get_22(), mass_dimension::minus3);

    return *this;
  }


std::string LSSEFT::make_unique_kernel_name()
  {
    return this->kernel_root + std::to_string(this->kernel_count++);
  }


void LSSEFT::process_kernels(const Pk_rsd_group& group, LSSEFT_impl::mass_dimension dim)
  {
    auto visitor = [&](const one_loop_element& elt) -> void
      {
        using LSSEFT_impl::LSSEFT_kernel;

        LSSEFT_kernel ker{elt.get_integrand(), elt.get_measure(), elt.get_Wick_product(),
                          elt.get_integration_variables(), elt.get_external_momenta(), dim};

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


    // KERNELS

    // generate container class
    this->write_container_class();

    // generate 'missing' statements for kernels
    this->write_kernel_missing();

    // generate store statements for kernels
    this->write_kernel_store();

    // generate find statements for kernels
    this->write_kernel_find();

    // generate kernel integrands
    this->write_kernel_integrands();

    // write kernel integrate statements
    this->write_integrate_stmts();

    // write kernel drop-index statements
    this->write_kernel_dropidx_stmts();

    // write kernel make-index statements
    this->write_kernel_makeidx_stmts();


    // ONE LOOP POWER SPECTRA

    // generate 'missing' statements for one-loop Pks
    this->write_Pk_missing();

    // generate store statements for kernels
    this->write_Pk_store();

    // write compute statements for the different Pk
    this->write_Pk_compute_stmts();

    // write find statements for the different Pk
    this->write_Pk_find();

    // write compute functions for the different Pk
    this->write_Pk_expressions();

    //! write Pk drop-index statements
    this->write_Pk_dropidx_stmts();

    //! write Pk make-index statements
    this->write_Pk_makeidx_stmts();


    // MULTIPOLE POWER SPECTRA

    // generate 'missing' statements for one-loop Pks
    this->write_multipole_missing();

    // write 'decompose' statements to construct Legendre modes
    this->write_multipole_decompose_stmts();

    // write store statements for Pn
    this->write_multipole_store();

    //! write multipole drop-index statements
    this->write_multipole_dropidx_stmts();

    //! write multiple mkake-index statements
    this->write_multipole_makeidx_stmts();
  }


boost::filesystem::path LSSEFT::make_output_path(const boost::filesystem::path& leaf) const
  {
    auto path = this->root / boost::filesystem::path{"autogenerated"} / leaf;
    return path;
  }


void LSSEFT::write_create() const
  {
    auto output = this->make_output_path("create_stmts.cpp");

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


void LSSEFT::write_kernel_integrands() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("kernel_integrands.cpp");

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

        const auto& integration_vars = kernel.get_integration_variables();
        const auto& external_momenta = kernel.get_external_momenta();
        const auto& k = *external_momenta.begin();

        GiNaC::exmap subs_map = { {q0, q_}, {x, z_}, {k, k_} };

        // write create statements for all kernels that we require
        outf << "static int " << name << "_integrand(const int* ndim_, const cubareal x_[], const int* ncomp_, cubareal f_[], void* userdata_)" << '\n';
        outf << " {" << '\n';
        outf << "   using oneloop_momentum_impl::integrand_data;" << '\n';
        outf << "   integrand_data* data_ = static_cast<integrand_data*>(userdata_);" << '\n';
        outf << '\n';
        outf << "   double k_ = data_->k * Mpc_units::Mpc;" << '\n';
        outf << "   double q_ = (data_->IR_cutoff + x_[0] * data_->q_range) * Mpc_units::Mpc;" << '\n';

        bool has_x_integral = false;
        if(integration_vars.find(x) != integration_vars.end())
          {
            has_x_integral = true;
            outf << "   double z_ = 2.0*x_[1] - 1.0;" << '\n';
          }
        else
          {
            outf << "   // no z_ integral in this kernel; measure should be jacobian_dq" << '\n';
          }

        outf << '\n';
        outf << "   double value_ = " << kernel.print_integrand(subs_map) << ";" << '\n';
        outf << "   double measure_ = " << kernel.print_measure(subs_map) << ";" << '\n';
        outf << "   double Wick_ = " << kernel.print_WickProduct(subs_map, external_momenta) << ";" << '\n';
        outf << "   f_[0] = (";

        if(has_x_integral)
          {
            outf << "data_->jacobian_dqdx";
          }
        else
          {
            outf << "data_->jacobian_dq";
          }

        outf << " * Mpc_units::Mpc) * value_ * measure_ * Wick_;" << '\n';
        outf << '\n';
        outf << "   return 0;" << '\n';
        outf << " }" << '\n';
        outf << '\n';
      }

    outf.close();
  }


const std::map< LSSEFT_impl::mass_dimension, std::string > integral_type_map
  = { { LSSEFT_impl::mass_dimension::zero, "dimless_integral" },
      { LSSEFT_impl::mass_dimension::minus3, "inverse_energy3_integral" } };


const std::map< LSSEFT_impl::mass_dimension, std::string > integral_1322_map
  = { { LSSEFT_impl::mass_dimension::zero, "loop_integral_type::P22" },
      { LSSEFT_impl::mass_dimension::minus3, "loop_integral_type::P13" } };


void LSSEFT::write_container_class() const
  {
    using LSSEFT_impl::LSSEFT_kernel;
    using LSSEFT_impl::mass_dimension;

    auto output = this->make_output_path("kernel_class.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    outf << "class kernels" << '\n';
    outf << " {" << '\n';
    outf << '\n';

    outf << " public:" << '\n';

    // value constructor
    outf << '\n';
    outf << "    //! value constructor" << '\n';
    outf << "    kernels(";

    // constructor argument list
    unsigned int count = 0;
    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        mass_dimension dim = kernel.get_dimension();

        if(count > 0) outf << ", ";
        outf << "const " << integral_type_map.at(dim) << "& " << name << "_";
        ++count;
      }
    outf << ")" << '\n';

    // constructor initializer list
    outf << "     : fail(false)";
    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        outf << ", " << name << "(" << name << "_)";
      }
    outf << '\n';
    outf << "     {" << '\n';
    outf << "     }" << '\n';

    // empty constructor
    outf << '\n';
    outf << "    //! empty constructor" << '\n';
    outf << "    kernels()" << '\n';
    outf << "     : fail(false)";
    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        outf << ", " << name << "()";
      }
    outf << '\n';
    outf << "     {" << '\n';
    outf << "     }" << '\n';

    // destructor
    outf << '\n';
    outf << "   //! destructor is default" << '\n';
    outf << "   ~kernels() = default;" << '\n';

    outf << '\n';
    outf << "public:" << '\n';

    outf << "    //! get failure state" << '\n'
         << "    bool get_fail() const { return this->fail; }" << '\n'
         << '\n'
         << "    //! set failed flag" << '\n'
         << "    void mark_failed() { this->fail = true; }" << '\n';

    // accessors
    outf << '\n';
    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        mass_dimension dim = kernel.get_dimension();

        outf << "    " << integral_type_map.at(dim) << "& get_" << name << "() { return this->" << name << "; }" << '\n';
        outf << "    const " << integral_type_map.at(dim) << "& get_" << name << "() const { return this->" << name << "; }" << '\n';
      }

    // internal data
    outf << '\n';
    outf << " private:" << '\n';

    outf << '\n';
    outf << "    bool fail;" << '\n';

    outf << '\n';
    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        mass_dimension dim = kernel.get_dimension();

        outf << "    " << integral_type_map.at(dim) << " " << name << ";" << '\n';
      }

    outf << '\n';
    outf << "    friend class boost::serialization::access;" << '\n';

    outf << '\n';
    outf << "    template <typename Archive>" << '\n'
         << "    void serialize(Archive& ar, unsigned int version)" << '\n'
         << "     {" << '\n'
         << "       ar & fail;" << '\n';
    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        outf << "       ar & " << name << ";" << '\n';
      }
    outf << "     }" << '\n';

    // close class brace
    outf << " };" << '\n';
  }


void LSSEFT::write_integrate_stmts() const
  {
    using LSSEFT_impl::mass_dimension;
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("integrate_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    outf << "    kernels ker;" << '\n';
    outf << "    bool fail = false;" << '\n';

    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        mass_dimension dim = kernel.get_dimension();

        outf << "    fail |= this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::" << name
             << "_integrand, ker.get_" << name << "(), " << integral_1322_map.at(dim) << ", \"" << name << "\");" << '\n';
      }

    outf << '\n';
    outf << "    if(fail) ker.mark_failed();" << '\n';

    outf.close();
  }


void LSSEFT::write_kernel_store() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("store_kernel_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        outf << "store_impl::store_loop_kernel(db, \"" << name << "\", ker.get_" << name << "(), model, params, sample);" << '\n';
      }

    outf.close();
  }


void LSSEFT::write_kernel_missing() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("missing_kernel_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        outf << "loop_configs " << name
             << " = update_missing_loop_integral_configurations(db, model, params, Pk_lin, \"" << name
             << "\", required_configs, total_missing);" << '\n';
      }
    outf << '\n';

    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        outf << "drop_inconsistent_configurations(db, model, params, Pk_lin, \"" << name << "\", " << name << ", total_missing);" << '\n';
      }

    outf.close();
  }


void LSSEFT::write_kernel_find() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("find_kernel_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    outf << "kernels ker;" << '\n';

    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        outf << "find_impl::read_loop_kernel(db, \"" << name << "\", model, params, k, Pk, UV_cutoff, ker.get_" << name << "(), IR_cutoff);" << '\n';

      }

    outf.close();
  }


void LSSEFT::write_Pk_missing() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("missing_Pk_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        outf << "std::set<unsigned int> " << name
             << "_mu0 = update_missing_one_loop_Pk(db, model, growth_params, loop_params, init_Pk, final_Pk, \"" << name
             << "_mu0\", z_table, record, missing);" << '\n';

        outf << "std::set<unsigned int> " << name
             << "_mu2 = update_missing_one_loop_Pk(db, model, growth_params, loop_params, init_Pk, final_Pk, \"" << name
             << "_mu2\", z_table, record, missing);" << '\n';

        outf << "std::set<unsigned int> " << name
             << "_mu4 = update_missing_one_loop_Pk(db, model, growth_params, loop_params, init_Pk, final_Pk, \"" << name
             << "_mu4\", z_table, record, missing);" << '\n';

        outf << "std::set<unsigned int> " << name
             << "_mu6 = update_missing_one_loop_Pk(db, model, growth_params, loop_params, init_Pk, final_Pk, \"" << name
             << "_mu6\", z_table, record, missing);" << '\n';

        outf << "std::set<unsigned int> " << name
             << "_mu8 = update_missing_one_loop_Pk(db, model, growth_params, loop_params, init_Pk, final_Pk, \"" << name
             << "_mu8\", z_table, record, missing);" << '\n';
      }
    outf << '\n';

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        outf << "drop_inconsistent_redshifts(db, model, init_Pk, growth_params, loop_params, final_Pk, \"" << name
             << "_mu0\", record, " << name << "_mu0, missing);" << '\n';

        outf << "drop_inconsistent_redshifts(db, model, init_Pk, growth_params, loop_params, final_Pk, \"" << name
             << "_mu2\", record, " << name << "_mu2, missing);" << '\n';

        outf << "drop_inconsistent_redshifts(db, model, init_Pk, growth_params, loop_params, final_Pk, \"" << name
             << "_mu4\", record, " << name << "_mu4, missing);" << '\n';

        outf << "drop_inconsistent_redshifts(db, model, init_Pk, growth_params, loop_params, final_Pk, \"" << name
             << "_mu6\", record, " << name << "_mu6, missing);" << '\n';

        outf << "drop_inconsistent_redshifts(db, model, init_Pk, growth_params, loop_params, final_Pk, \"" << name
             << "_mu8\", record, " << name << "_mu8, missing);" << '\n';
      }

    outf.close();
  }


void LSSEFT::write_Pk_compute_stmts() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("compute_Pk_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        outf << "rsd_dd_Pk " << name << "_mu0 = compute_" << name << "_mu0(k, val.second, loop_data, Ptr_init, Ptr_final);" << '\n';
        outf << "rsd_dd_Pk " << name << "_mu2 = compute_" << name << "_mu2(k, val.second, loop_data, Ptr_init, Ptr_final);" << '\n';
        outf << "rsd_dd_Pk " << name << "_mu4 = compute_" << name << "_mu4(k, val.second, loop_data, Ptr_init, Ptr_final);" << '\n';
        outf << "rsd_dd_Pk " << name << "_mu6 = compute_" << name << "_mu6(k, val.second, loop_data, Ptr_init, Ptr_final);" << '\n';
        outf << "rsd_dd_Pk " << name << "_mu8 = compute_" << name << "_mu8(k, val.second, loop_data, Ptr_init, Ptr_final);" << '\n';

        outf << "Pks.emplace(std::make_pair(\"" << name
             << "\", oneloop_Pk{k_tok, gf_factors.get_params_token(), loop_data.get_params_token(), Pk_init.get_token(),final_tok, loop_data.get_IR_token(), loop_data.get_UV_token(), val.first.get_id(), "
             << name << "_mu0, " << name << "_mu2, " << name << "_mu4, " << name << "_mu6, " << name << "_mu8}));"
             << '\n';
        outf << '\n';
      }

    outf.close();
  }


void LSSEFT::write_Pk_store() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("store_Pk_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        outf << "const oneloop_Pk& " << name << " = record.at(\"" << name << "\");" << '\n';
        outf << "store_impl::store_one_loop_rsd_Pk(db, \"" << name << "_mu0\", " << name << ".get_dd_rsd_mu0(), model, " << name << ");" << '\n';
        outf << "store_impl::store_one_loop_rsd_Pk(db, \"" << name << "_mu2\", " << name << ".get_dd_rsd_mu2(), model, " << name << ");" << '\n';
        outf << "store_impl::store_one_loop_rsd_Pk(db, \"" << name << "_mu4\", " << name << ".get_dd_rsd_mu4(), model, " << name << ");" << '\n';
        outf << "store_impl::store_one_loop_rsd_Pk(db, \"" << name << "_mu6\", " << name << ".get_dd_rsd_mu6(), model, " << name << ");" << '\n';
        outf << "store_impl::store_one_loop_rsd_Pk(db, \"" << name << "_mu8\", " << name << ".get_dd_rsd_mu8(), model, " << name << ");" << '\n';
        outf << '\n';
      }

    outf.close();
  }


void LSSEFT::write_Pk_find() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("find_Pk_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        outf << "rsd_dd_Pk " << name << "_mu0;" << '\n';
        outf << "rsd_dd_Pk " << name << "_mu2;" << '\n';
        outf << "rsd_dd_Pk " << name << "_mu4;" << '\n';
        outf << "rsd_dd_Pk " << name << "_mu6;" << '\n';
        outf << "rsd_dd_Pk " << name << "_mu8;" << '\n';

        outf << "find_impl::read_dd_rsd_Pk(db, \"" << name << "_mu0\", model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, " << name << "_mu0);" << '\n';
        outf << "find_impl::read_dd_rsd_Pk(db, \"" << name << "_mu2\", model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, " << name << "_mu2);" << '\n';
        outf << "find_impl::read_dd_rsd_Pk(db, \"" << name << "_mu4\", model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, " << name << "_mu4);" << '\n';
        outf << "find_impl::read_dd_rsd_Pk(db, \"" << name << "_mu6\", model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, " << name << "_mu6);" << '\n';
        outf << "find_impl::read_dd_rsd_Pk(db, \"" << name << "_mu8\", model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, " << name << "_mu8);" << '\n';

        outf << "payload->emplace(std::make_pair(\"" << name
             << "\", oneloop_Pk{k, growth_params, loop_params, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, z, "
             << name << "_mu0, " << name << "_mu2, " << name << "_mu4, " << name << "_mu6, " << name << "_mu8}));"
             << '\n';

        outf << '\n';
      }
  }


void LSSEFT::write_Pk_mu_component(std::ofstream& outf, const std::string& name, const Pk_rsd& Pk, unsigned int mu) const
  {
    std::string tag = std::string{"mu"} + std::to_string(mu);

    outf << "rsd_dd_Pk compute_" << name << "_" << tag
         << "(const Mpc_units::energy& k, const oneloop_growth_record& val, const loop_integral& loop_data, const Pk_value& Ptr_init, const boost::optional<Pk_value>& Ptr_final)"
         << '\n';

    outf << " {" << '\n';

    outf << "   const kernels& ker = loop_data.get_kernels();" << '\n';
    outf << '\n';

    const auto& tree = Pk.get_tree();
    outf << "   Pk_value tree";
    std::ostringstream tree_buffer;
    unsigned int count = 0;
    auto tree_writer = [&](const one_loop_element& elt) -> void
      {
        using LSSEFT_impl::format_print;

        if(count > 0) tree_buffer << " + ";
        tree_buffer << format_print(elt.get_integrand()*elt.get_time_function());    // have to include "integrand" at tree-level, which is really a normalization factor
        ++count;
      };
    tree.visit({mu}, tree_writer);
    if(count == 0) outf << ";   // no contribution at mu^" << mu;
    else           outf << " = (" << tree_buffer.str() << ") * (Ptr_final ? *Ptr_final : Ptr_init);";
    outf << '\n' << '\n';

    const auto& P13 = Pk.get_13();
    outf << "   Pk_value P13";
    std::ostringstream P13_buffer;
    count = 0;
    auto P13_writer = [&](const one_loop_element& elt) -> void
      {
        using LSSEFT_impl::LSSEFT_kernel;
        using LSSEFT_impl::mass_dimension;
        using LSSEFT_impl::format_print;

        LSSEFT_kernel ker{elt.get_integrand(), elt.get_measure(), elt.get_Wick_product(),
                          elt.get_integration_variables(), elt.get_external_momenta(), mass_dimension::zero};

        const std::string& kname = this->kernel_db.at(ker);

        if(count > 0) P13_buffer << " + ";
        P13_buffer << "(" << format_print(elt.get_time_function()) << ")*ker.get_" << kname << "()";
        ++count;
      };
    P13.visit({mu}, P13_writer);
    if(count == 0) outf << ";   // no contribution at mu^" << mu;
    else           outf << " = Ptr_init * (" << P13_buffer.str() << ");";
    outf << '\n' << '\n';

    const auto& P22 = Pk.get_22();
    outf << "   Pk_value P22";
    std::ostringstream P22_buffer;
    count = 0;
    auto P22_writer = [&](const one_loop_element& elt) -> void
      {
        using LSSEFT_impl::LSSEFT_kernel;
        using LSSEFT_impl::mass_dimension;
        using LSSEFT_impl::format_print;

        LSSEFT_kernel ker{elt.get_integrand(), elt.get_measure(), elt.get_Wick_product(),
                          elt.get_integration_variables(), elt.get_external_momenta(), mass_dimension::minus3};

        const std::string& kname = this->kernel_db.at(ker);

        if(count > 0) P22_buffer << " + ";
        P22_buffer << "(" << format_print(elt.get_time_function()) << ")*ker.get_" << kname << "()";
        ++count;
      };
    P22.visit({mu}, P22_writer);
    if(count == 0) outf << ";   // no contribution at mu^" << mu;
    else           outf << " = " << P22_buffer.str() << ";";
    outf << '\n' << '\n';

    outf << "   return rsd_dd_Pk{tree, P13, P22};" << '\n';

    outf << " }" << '\n';
    outf << '\n';
  }


void LSSEFT::write_Pk_expressions() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("Pk_expressions.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;
        const Pk_rsd& Pk = record.second;

        this->write_Pk_mu_component(outf, name, Pk, 0);
        this->write_Pk_mu_component(outf, name, Pk, 2);
        this->write_Pk_mu_component(outf, name, Pk, 4);
        this->write_Pk_mu_component(outf, name, Pk, 6);
        this->write_Pk_mu_component(outf, name, Pk, 8);
      }

    outf.close();
  }


void LSSEFT::write_kernel_dropidx_stmts() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("dropidx_kernel_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        outf << "sqlite3_operations::drop_index(this->handle, \"" << name << "\", { \"mid\", \"params_id\", \"kid\", \"Pk_id\", \"IR_id\", \"UV_id\" });" << '\n';
      }

    outf.close();
  }


void LSSEFT::write_kernel_makeidx_stmts() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("makeidx_kernel_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->kernel_db)
      {
        const LSSEFT_kernel& kernel = record.first;
        const std::string& name = record.second;

        outf << "sqlite3_operations::create_index(this->handle, \"" << name << "\", { \"mid\", \"params_id\", \"kid\", \"Pk_id\", \"IR_id\", \"UV_id\" });" << '\n';
      }

    outf.close();
  }


void LSSEFT::write_Pk_dropidx_stmts() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("dropidx_Pk_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;
        const Pk_rsd& Pk = record.second;

        outf << "sqlite3_operations::drop_index(this->handle, \"" << name
             << "_mu0\", { \"mid\", \"growth_params\", \"loop_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_id\", \"UV_id\", });"
             << '\n';

        outf << "sqlite3_operations::drop_index(this->handle, \"" << name
             << "_mu2\", { \"mid\", \"growth_params\", \"loop_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_id\", \"UV_id\", });"
             << '\n';

        outf << "sqlite3_operations::drop_index(this->handle, \"" << name
             << "_mu4\", { \"mid\", \"growth_params\", \"loop_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_id\", \"UV_id\", });"
             << '\n';

        outf << "sqlite3_operations::drop_index(this->handle, \"" << name
             << "_mu6\", { \"mid\", \"growth_params\", \"loop_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_id\", \"UV_id\", });"
             << '\n';

        outf << "sqlite3_operations::drop_index(this->handle, \"" << name
             << "_mu8\", { \"mid\", \"growth_params\", \"loop_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_id\", \"UV_id\", });"
             << '\n';
      }

    outf.close();
  }


void LSSEFT::write_Pk_makeidx_stmts() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("makeidx_Pk_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        outf << "sqlite3_operations::create_index(this->handle, \"" << name
             << "_mu0\", { \"mid\", \"growth_params\", \"loop_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_id\", \"UV_id\", });"
             << '\n';

        outf << "sqlite3_operations::create_index(this->handle, \"" << name
             << "_mu2\", { \"mid\", \"growth_params\", \"loop_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_id\", \"UV_id\", });"
             << '\n';

        outf << "sqlite3_operations::create_index(this->handle, \"" << name
             << "_mu4\", { \"mid\", \"growth_params\", \"loop_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_id\", \"UV_id\", });"
             << '\n';

        outf << "sqlite3_operations::create_index(this->handle, \"" << name
             << "_mu6\", { \"mid\", \"growth_params\", \"loop_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_id\", \"UV_id\", });"
             << '\n';

        outf << "sqlite3_operations::create_index(this->handle, \"" << name
             << "_mu8\", { \"mid\", \"growth_params\", \"loop_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_id\", \"UV_id\", });"
             << '\n';

        outf << '\n';
      }

    outf.close();
  }


void LSSEFT::write_multipole_missing() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("missing_multipole_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        outf << "std::set<unsigned int> " << name
             << "_P0 = update_missing_multipole_Pk(db, model, growth_params, loop_params, XY_params, init_Pk, final_Pk, \""
             << name << "_P0\", z_table, record, missing);" << '\n';

        outf << "std::set<unsigned int> " << name
             << "_P2 = update_missing_multipole_Pk(db, model, growth_params, loop_params, XY_params, init_Pk, final_Pk, \""
             << name << "_P2\", z_table, record, missing);" << '\n';

        outf << "std::set<unsigned int> " << name
             << "_P4 = update_missing_multipole_Pk(db, model, growth_params, loop_params, XY_params, init_Pk, final_Pk, \""
             << name << "_P4\", z_table, record, missing);" << '\n';
      }
    outf << '\n';

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        outf << "drop_inconsistent_redshifts(db, model, growth_params, loop_params, XY_params, init_Pk, final_Pk, \""
             << name << "_P0\", record, " << name << "_P0, missing);" << '\n';

        outf << "drop_inconsistent_redshifts(db, model, growth_params, loop_params, XY_params, init_Pk, final_Pk, \""
             << name << "_P2\", record, " << name << "_P2, missing);" << '\n';

        outf << "drop_inconsistent_redshifts(db, model, growth_params, loop_params, XY_params, init_Pk, final_Pk, \""
             << name << "_P4\", record, " << name << "_P4, missing);" << '\n';
      }

    outf.close();
  }


void LSSEFT::write_multipole_decompose_stmts() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("multipole_compute_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        outf << "apply(\"" << name << "\");" << '\n';
      }

    outf.close();
  }


void LSSEFT::write_multipole_store() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("store_multipole_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;

        outf << "const multipole_Pk& " << name << " = sample.at(\"" << name << "\");" << '\n';
        outf << "store_impl::store_multipole_Pk(db, \"" << name << "_P0\", " << name << ".get_P0(), model, " << name << ");" << '\n';
        outf << "store_impl::store_multipole_Pk(db, \"" << name << "_P2\", " << name << ".get_P2(), model, " << name << ");" << '\n';
        outf << "store_impl::store_multipole_Pk(db, \"" << name << "_P4\", " << name << ".get_P4(), model, " << name << ");" << '\n';
        outf << '\n';
      }

    outf.close();
  }


void LSSEFT::write_multipole_dropidx_stmts() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("dropidx_multipole_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;
        const Pk_rsd& Pk = record.second;

        outf << "sqlite3_operations::drop_index(this->handle, \"" << name
             << "_P0\", { \"mid\", \"growth_params\", \"loop_params\", \"XY_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_cutoff_id\", \"UV_cutoff_id\", \"IR_resum_id\" });"
             << '\n';

        outf << "sqlite3_operations::drop_index(this->handle, \"" << name
             << "_P2\", { \"mid\", \"growth_params\", \"loop_params\", \"XY_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_cutoff_id\", \"UV_cutoff_id\", \"IR_resum_id\" });"
             << '\n';

        outf << "sqlite3_operations::drop_index(this->handle, \"" << name
             << "_P4\", { \"mid\", \"growth_params\", \"loop_params\", \"XY_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_cutoff_id\", \"UV_cutoff_id\", \"IR_resum_id\" });"
             << '\n';
      }

    outf.close();
  }


void LSSEFT::write_multipole_makeidx_stmts() const
  {
    using LSSEFT_impl::LSSEFT_kernel;

    auto output = this->make_output_path("makeidx_multipole_stmts.cpp");

    std::ofstream outf{output.string(), std::ios_base::out | std::ios_base::trunc};

    for(const auto& record : this->Pk_db)
      {
        const std::string& name = record.first;
        const Pk_rsd& Pk = record.second;

        outf << "sqlite3_operations::create_index(this->handle, \"" << name
             << "_P0\", { \"mid\", \"growth_params\", \"loop_params\", \"XY_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_cutoff_id\", \"UV_cutoff_id\", \"IR_resum_id\" });"
             << '\n';

        outf << "sqlite3_operations::create_index(this->handle, \"" << name
             << "_P2\", { \"mid\", \"growth_params\", \"loop_params\", \"XY_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_cutoff_id\", \"UV_cutoff_id\", \"IR_resum_id\" });"
             << '\n';

        outf << "sqlite3_operations::create_index(this->handle, \"" << name
             << "_P4\", { \"mid\", \"growth_params\", \"loop_params\", \"XY_params\", \"kid\", \"zid\", \"init_Pk_id\", \"final_Pk_id\", \"IR_cutoff_id\", \"UV_cutoff_id\", \"IR_resum_id\" });"
             << '\n';
      }

    outf.close();
  }
