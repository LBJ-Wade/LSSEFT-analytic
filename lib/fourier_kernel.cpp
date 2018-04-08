//
// Created by David Seery on 21/08/2017.
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

#include <algorithm>
#include <sstream>
#include <string>

#include "fourier_kernel.h"
#include "lib/correlators/detail/relabel_product.h"
#include "lib/correlators/detail/Rayleigh_momenta.h"
#include "utilities/GiNaC_utils.h"


namespace fourier_kernel_impl
  {

    key::key(kernel& k)
      : tm(k.tm),
        iv(k.iv)
      {
      }


    key::key(const time_function& tm_, const initial_value_set& iv_)
      : tm(tm_),
        iv(iv_)
      {
      }


    std::vector<GiNaC::symbol>
    key::get_ordered_iv_symbols() const
      {
        std::vector<GiNaC::symbol> symbols;
        for(auto t = this->iv.value_cbegin(); t != this->iv.value_cend(); ++t)
          {
            symbols.push_back(t->get_symbol());
          }
        
        std::sort(symbols.begin(), symbols.end(), std::less<GiNaC::symbol>{});
        
        return symbols;
      }


    size_t key::hash() const
      {
        // to hash the time expression, expand it completely and print
        // we are guaranteed that the expressions comes in a canonical order, even if that order is unpredictable
        std::ostringstream time_string;
        time_string << this->tm.expand();

        size_t h = 0;
        hash_impl::hash_combine(h, time_string.str());

        // to hash the initial value set, order its symbols lexicographically
        const auto symbols = this->get_ordered_iv_symbols();

        std::string symbol_string;
        std::for_each(symbols.begin(), symbols.end(),
                      [&](const GiNaC::symbol& e) -> void
                        { symbol_string += e.get_name(); });

        // combine both hashes together
        hash_impl::hash_combine(h, symbol_string);

        // return final value
        return h;
      }


    bool key::is_equal(const key& obj) const
      {
        const auto& at = this->tm;
        const auto& bt = obj.tm;

        // test for equality of expressions
        auto rt = (at == bt);
        if(!static_cast<bool>(rt)) return false;

        // test for equality of initial-value strings
        // we do this by ordering their symbol names lexicographically
        // and testing for equality of those
        auto a_symbols = this->get_ordered_iv_symbols();
        auto b_symbols = obj.get_ordered_iv_symbols();

        return std::equal(a_symbols.cbegin(), a_symbols.cend(),
                          b_symbols.cbegin(), b_symbols.cend(),
                          [](const GiNaC::symbol& asym, const GiNaC::symbol& bsym) -> bool
                            { return asym.get_name() == bsym.get_name(); });
      }


    kernel::kernel(GiNaC::ex K_, initial_value_set iv_, time_function tm_, subs_list vs_, service_locator& lc_)
      : K(std::move(K_)),
        tm(std::move(tm_)),
        iv(std::move(iv_)),
        vs(std::move(vs_)),
        loc(lc_)
      {
        // normalize the time function, redistributing factors into the kernel if needed
        auto norm = get_normalization_factor(tm, loc);
        tm /= norm;
        K *= norm;
      }


    kernel::kernel(initial_value_set iv_, service_locator& sl_)
      : K(1),
        tm(1),
        iv(std::move(iv_)),
        loc(sl_)
      {
      }

    kernel& kernel::operator+=(const kernel& rhs)
      {
        // kernels can only be added if their time functions agree
        if(static_cast<bool>(this->tm != rhs.tm))
          throw exception(ERROR_CANNOT_ADD_KERNELS_WITH_UNEQUAL_TIME_FUNCTIONS, exception_code::loop_transformation_error);

        // we need to relabel momenta in the right-hand-side if they do not match
        // our current labelling
        auto our_mma = this->get_ordered_momenta();
        auto their_mma = rhs.get_ordered_momenta();
        
        // check that momenta are compatible in the sense that their (ordered) symbol names agree
        if(!std::equal(our_mma.cbegin(), our_mma.cend(),
                       their_mma.cbegin(), their_mma.cend(),
                       [](const auto& a, const auto& b) -> bool
                         { return a.get().get_symbol().get_name() == b.get().get_symbol().get_name(); }))
          throw exception(ERROR_KERNEL_INITIAL_VALUES_DISAGREE, exception_code::kernel_error);
        
        // build a substitution map for those momenta that disagree
        GiNaC::exmap mma_map;
        
        for(auto u = our_mma.cbegin(), v = their_mma.cbegin(); u != our_mma.cend() && v != their_mma.cend(); ++u, ++v)
          {
            const auto& our_q = u->get().get_momentum();
            const auto& their_q = v->get().get_momentum();

            // if the momentum labels disagree, then remap theirs
            if(our_q != their_q) mma_map[their_q] = our_q;
          }
        
        // also need to remap and merge any substitution rules in the right-hand side
        using Rayleigh::merge_Rayleigh_rules;
        auto relabel_map = merge_Rayleigh_rules(this->vs, rhs.vs, this->iv.get_momenta(), mma_map, this->loc);
        std::copy(relabel_map.begin(), relabel_map.end(), std::inserter(mma_map, mma_map.begin()));

        // now perform relabelling in the kernel
        // note there's no need to perform *index* relabelling since this is a sum -- relabelling indices
        // is needed only in a product
        auto temp = simplify_index(this->K + rhs.K.subs(mma_map), this->vs, this->loc);
        this->K = temp;

        return *this;
      }
    
    
    kernel& kernel::operator*=(const kernel& rhs)
      {
        auto& sf = this->loc.get_symbol_factory();

        // product time function is just product of each individual time function
        this->tm *= rhs.tm;
    
        // build a substitution rule for all momenta in rhs that also occur for us
        GiNaC::exmap mma_map;
        
        const auto our_syms = this->iv.get_momenta();
        for(auto t = rhs.iv.value_cbegin(); t != rhs.iv.value_cend(); ++t)
          {
            const GiNaC::symbol& label = t->get_momentum();
            bool collision = false;
            if(our_syms.find(label) != our_syms.end()) collision = true;
            if(this->vs.find(label) != this->vs.end()) collision = true;
            
            // if no collision can just keep the old symbol
            if(!collision) continue;
            
            // otherwise, manufacture a new symbol
            auto relabel = sf.make_unique_momentum();
            
            mma_map[label] = relabel;
          }
    
        // merge RHS initial value list with ours
        for(auto t = rhs.iv.value_cbegin(); t != rhs.iv.value_cend(); ++t)
          {
            // get momentum for this field
            const auto& q = t->get_momentum();
            
            // is this momentum being relabelled?
            auto u = mma_map.find(q);
            if(u != mma_map.end())
              {
                // yes, so merge the relabelled version
                const auto& u_sym = GiNaC::ex_to<GiNaC::symbol>(u->second);
                this->iv.insert(t->relabel_momentum(u_sym));
              }
            else
              {
                // no, so merge the original
                this->iv.insert(*t);
              }
          }
    
        // merge RHS substitution list with ours
        // (occurs after merging initial value list so all reserved symbols are captured)
        using Rayleigh::merge_Rayleigh_rules;
        auto relabel_map = merge_Rayleigh_rules(this->vs, rhs.vs, this->iv.get_momenta(), mma_map, this->loc);
        std::copy(relabel_map.begin(), relabel_map.end(), std::inserter(mma_map, mma_map.begin()));

        // build final expression, performing any necessary index (or other) relabellings on RHS
        auto temp = simplify_index(relabel_index_product(this->K, rhs.K.subs(mma_map), this->loc), this->vs, this->loc);
        this->K = temp;

        // renormalize time function
        auto norm = get_normalization_factor(tm, loc);
        this->tm /= norm;
        this->K *= norm;

        return *this;
      }
    
    
    kernel::momenta_list kernel::get_ordered_momenta() const
      {
        momenta_list list;
        
        for(auto t = this->iv.value_cbegin(); t != this->iv.value_cend(); ++t)
          {
            list.emplace_back(std::cref(*t));
          }
    
        std::sort(list.begin(), list.end(),
                  [](const momenta_list::value_type& a, const momenta_list::value_type& b) -> bool
                    { return std::less<GiNaC::symbol>{}(a.get().get_symbol(), b.get().get_symbol()); });
        
        return list;
      }
    
    
    kernel operator-(const kernel& a)
      {
        kernel b = a;
        b.K = -b.K;
        return b;
      }
    
    
    kernel operator+(const kernel& a, const kernel& b)
      {
        kernel c = a;
        c += b;
        return c;
      }
    
    
    kernel operator-(const kernel& a, const kernel& b)
      {
        return a + (-b);
      }
    
    
    kernel operator*(const GiNaC::ex a, const kernel& b)
      {
        kernel c = b;

        // multiply kernel b by the GiNaC expression a
        // depending what is in a, we want to handle this operation differently
        // we'd like the timefunction to include only SPT time-dependent function and possibly numeric constants
        // we'd like everything else, such as parameters (eg. halo bias) to go into the integrand

        // explicitly distribute a over a top-level sum if one is present
        // this is because GiNaC is very reluctant to simplify cancelling factors in such cases,
        // meaning that we may get uncancelled factors of (eg.) 1+z in the result, even though they all cancel down to 1
        GiNaC::ex time_factor;
        GiNaC::ex integrand_factor;

        std::tie(time_factor, integrand_factor) = partition_factor(a, b.loc);

        c.tm = GiNaC::collect_common_factors(c.tm.expand() * time_factor);
        c.K *= integrand_factor;
        
        // adjust normalization of time function if needed
        auto norm = get_normalization_factor(c.tm, c.loc);
        c.tm /= norm;
        c.K *= norm;

        return c;
      }
    
    
    kernel operator*(const kernel& a, const GiNaC::ex b)
      {
        return b*a;
      }
    
    
    kernel operator*(const kernel& a, const kernel& b)
      {
        kernel c = a;
        c *= b;
        return c;
      }
    
    
    kernel operator/(const kernel& a, const GiNaC::ex b)
      {
        return (GiNaC::numeric{1}/b) * a;
      }
    
    
    kernel diff_z(const kernel& a)
      {
        auto& sf = a.loc.get_symbol_factory();

        kernel b = a;
        const auto& z = sf.get_z();
        b.tm = GiNaC::diff(b.tm, z);

        // renormalize time function
        auto norm = get_normalization_factor(b.tm, b.loc);
        b.tm /= norm;
        b.K *= norm;

        return b;
      }
    
    
    vector kernel::get_total_momentum() const
      {
        // we're guaranteed that the initial value list is nonempty, so it's safe to dereference its first element
        vector sum = *this->iv.value_cbegin();
        
        for(auto t = ++this->iv.value_cbegin(); t != this->iv.value_cend(); ++t)
          {
            sum += *t;
          }
        
        return sum;
      }
    
    
    kernel& kernel::multiply_kernel(GiNaC::ex f)
      {
        auto our_syms = this->iv.get_momenta();
        const auto& params = this->loc.get_symbol_factory().get_parameters();
        std::copy(params.begin(), params.end(), std::inserter(our_syms, our_syms.begin()));
    
        // ensure that f only involves symbols in the initial value set or s, or declared as parameters
        auto expr_syms = get_expr_symbols(f);
    
        for(const auto& sym : expr_syms)
          {
            if(our_syms.find(sym) != our_syms.end()) continue;
        
            std::ostringstream msg;
            msg << ERROR_MULTIPLY_KERNEL_UNKNOWN_MOMENTUM << " '" << sym << "'";
            throw exception(msg.str(), exception_code::kernel_error);
          }

        this->K *= f;

        // renormalize time function
        auto norm = get_normalization_factor(tm, loc);
        this->tm /= norm;
        this->K *= norm;

        return *this;
      }
    
    
    kernel& kernel::multiply_kernel(GiNaC::ex f, GiNaC::symbol s, GiNaC::ex rule)
      {
        auto our_syms = this->iv.get_momenta();
        const auto& params = this->loc.get_symbol_factory().get_parameters();
        std::copy(params.begin(), params.end(), std::inserter(our_syms, our_syms.begin()));

        // ensure that f only involves symbols in the initial value set or s
        auto expr_syms = get_expr_symbols(f);

        for(const auto& sym : expr_syms)
          {
            if(our_syms.find(sym) != our_syms.end()) continue;
            if(sym.get_name() == s.get_name()) continue;
            
            std::ostringstream msg;
            msg << ERROR_MULTIPLY_KERNEL_UNKNOWN_MOMENTUM << " '" << sym << "'";
            throw exception(msg.str(), exception_code::kernel_error);
          }
        
        // ensure that rule only involves symbols in the initial value set
        auto rule_syms = get_expr_symbols(rule);
        for(const auto& sym : rule_syms)
          {
            if(our_syms.find(sym) != our_syms.end()) continue;
            if(sym.get_name() == s.get_name()) continue;
    
            std::ostringstream msg;
            msg << ERROR_MULTIPLY_KERNEL_UNKNOWN_MOMENTUM << " '" << sym << "'";
            throw exception(msg.str(), exception_code::kernel_error);
          }
        
        // ensure that s doesn't already exist in replacement rule set
        if(this->vs.find(s) != this->vs.end())
          {
            std::ostringstream msg;
            msg << ERROR_SUBSTITUTION_RULE_ALREADY_EXISTS << " '" << s << "'";
            throw exception(msg.str(), exception_code::kernel_error);
          }

        // check whether another rule with the same (or equivalent) RHS already exists
        auto t = this->vs.begin();
        for(; t != this->vs.end(); ++t)
          {
            if(static_cast<bool>(t->second == rule))
              {
                GiNaC::exmap map;
                map[s] = t->first;
                f = f.subs(map);
                break;
              }

            if(static_cast<bool>(t->second == -rule))
              {
                GiNaC::exmap map;
                map[s] = -t->first;
                f = f.subs(map);
                break;
              }
          }
        
        this->K *= f;

        // insert new rule if needed
        if(t == this->vs.end()) this->vs[s] = rule;

        // renormalize time function
        auto norm = get_normalization_factor(tm, loc);
        this->tm /= norm;
        this->K *= norm;

        return *this;
      }
    
    
    void kernel::write(std::ostream& out) const
      {
        out << "  " << MESSAGE_KERNEL_TIME_FUNCTION << " = " << this->tm << '\n';
    
        out << "  " << MESSAGE_KERNEL_IVSET << " = ";
        for(auto u = this->iv.value_cbegin(); u != this->iv.value_cend(); ++u)
          {
            out << " " << u->get_symbol() << "(" << u->get_momentum() << ")";
          }
        out << '\n';
    
        if(!this->vs.empty())
          {
            out << "  " << MESSAGE_KERNEL_RAYLEIGH_RULES << ":" << '\n';
            for(const auto& u : vs)
              {
                out << "    " << u.first << " -> " << u.second << '\n';
              }
          }
    
        out << "  " << MESSAGE_KERNEL_EXPRESSION << " = " << this->K << '\n';
      }


    void kernel::to_EdS()
      {
        auto z = this->loc.get_symbol_factory().get_z();

        GiNaC::exmap map =
          {
            { SPT::DA(z), 3*SPT::D(z)*SPT::D(z)/7 },
            { SPT::DB(z), 2*SPT::D(z)*SPT::D(z)/7 },
            { SPT::DD(z), 2*SPT::D(z)*SPT::D(z)*SPT::D(z)/21 },
            { SPT::DE(z), 4*SPT::D(z)*SPT::D(z)*SPT::D(z)/63 },
            { SPT::DF(z), SPT::D(z)*SPT::D(z)*SPT::D(z)/14 },
            { SPT::DG(z), SPT::D(z)*SPT::D(z)*SPT::D(z)/21 },
            { SPT::DJ(z), SPT::D(z)*SPT::D(z)*SPT::D(z)/9 },
            { SPT::fA(z), 2*SPT::f(z) },
            { SPT::fB(z), 2*SPT::f(z) },
            { SPT::fD(z), 3*SPT::f(z) },
            { SPT::fE(z), 3*SPT::f(z) },
            { SPT::fF(z), 3*SPT::f(z) },
            { SPT::fG(z), 3*SPT::f(z) },
            { SPT::fJ(z), 3*SPT::f(z) },
          };

        this->tm = this->tm.subs(map);
      }


    std::ostream& operator<<(std::ostream& str, const kernel& a)
      {
        a.write(str);
        return str;
      }
    
  }   // namespace fourier_kernel_impl


bool validate_ivset_nonempty(const initial_value_set& s, const GiNaC::ex& K, bool silent)
  {
    // ensure initial value set s is not empty
    if(s.empty())
      {
        if(!silent)
          {
            error_handler err;
            err.warn(WARNING_ORDER_ZERO_KERNEL);
            
            std::ostringstream msg;
            msg << MESSAGE_KERNEL_EXPRESSION << " = " << K;
            err.info(msg.str());
          }
        
        return false;
      }
    
    return true;
  }
  

void validate_subslist(const initial_value_set& s, const subs_list& vs)
  {
    // ensure that substitution list is a map from plain symbols to expressions, and
    // ensure that this list shares no common symbols with the initial value set
    auto s_mma = s.get_momenta();
    for(const auto& ele : vs)
      {
        const GiNaC::ex& label = ele.first;
        if(!GiNaC::is_exactly_a<GiNaC::symbol>(label))
          {
            std::ostringstream msg;
            msg << ERROR_SUBSTITION_LABEL_NOT_A_SYMBOL << " '" << label << "'";
            throw exception(msg.str(), exception_code::kernel_error);
          }
        
        const GiNaC::symbol& sym = GiNaC::ex_to<GiNaC::symbol>(label);
        if(s_mma.find(sym) != s_mma.end())
          {
            std::ostringstream msg;
            msg << ERROR_SUBSTITUTION_LIST_HAS_IV_MOMENTUM << " '" << sym << "'";
            throw exception(msg.str(), exception_code::kernel_error);
          }
      }
  }


void validate_structure(const GiNaC::ex& K)
  {
    auto idcs = K.get_free_indices();
    if(!idcs.empty())
      {
        throw exception(ERROR_KERNEL_NOT_SCALAR, exception_code::kernel_error);
      }
    
    // validate that the kernel K is a rational function
    // more general kernels are not supported (yet)
    if(!is_rational(K))
      {
        std::cerr << K << '\n';
        throw exception(ERROR_KERNEL_NOT_RATIONAL, exception_code::kernel_error);
      }
  }


void validate_momenta(const initial_value_set& s, const subs_list& vs, const GiNaC::ex& K,
                      const GiNaC_symbol_set& params, bool silent)
  {
    auto used = get_expr_symbols(K);
    auto avail = s.get_momenta();
    auto avail_plus_params = avail;
    
    // insert any symbols from the substitution list vs into the available set
    std::for_each(vs.begin(), vs.end(), [&](const subs_list::value_type& v) -> void
      {
        auto sym = GiNaC::ex_to<GiNaC::symbol>(v.first);
        avail_plus_params.insert(sym);
      });

    // also insert any symbols that have been declared as parameters
    std::copy(params.begin(), params.end(), std::inserter(avail_plus_params, avail_plus_params.begin()));
    
    GiNaC_symbol_set used_not_avail;
    GiNaC_symbol_set avail_not_used;
    
    // perform set differencing to find mismatch between available and used symbols
    std::set_difference(used.cbegin(), used.cend(), avail_plus_params.cbegin(), avail_plus_params.cend(),
                        std::inserter(used_not_avail, used_not_avail.end()), std::less<GiNaC::symbol>{});
    std::set_difference(avail.cbegin(), avail.cend(), used.cbegin(), used.cend(),
                        std::inserter(avail_not_used, avail_not_used.end()), std::less<GiNaC::symbol>{});
    
    // raise exception on attempt to use a Fourier momentum that isn't available
    if(!used_not_avail.empty())
      {
        std::cerr << "Kernel K = " << K << '\n';
        std::cerr << "Initial values =";
        for(auto t = s.value_cbegin(); t != s.value_cend(); ++t)
          {
            std::cerr << " " << t->get_symbol() << "(" << t->get_momentum() << ")";
          }
        std::cerr << '\n';
        
        std::ostringstream msg;
        msg << (used_not_avail.size() == 1 ? ERROR_UNKNOWN_MOMENTA_SING : ERROR_UNKNOWN_MOMENTA_PLURAL) << " ";
        
        // attach unknown momenta to error message
        unsigned int count = 0;
        for(const auto& k : used_not_avail)
          {
            if(count > 0) msg << ", ";
            msg << k;
            ++count;
          }
        
        // throw
        throw exception(msg.str(), exception_code::kernel_error);
      }
    
    // issue warning if there are available momenta that aren't used in the kernel
    if(!avail_not_used.empty() && !silent)
      {
        std::ostringstream msg;
        msg << (avail_not_used.size() == 1 ? WARNING_UNUSED_MOMENTA_SING : WARNING_UNUSED_MOMENTA_PLURAL) << " ";
        
        // attach unused momenta to message
        unsigned int count = 0;
        for(const auto& k : avail_not_used)
          {
            if(count > 0) msg << ", ";
            msg << k;
            ++count;
          }
        
        // issue error message
        error_handler err;
        err.warn(msg.str());
        
        std::ostringstream ker;
        ker << MESSAGE_KERNEL_EXPRESSION << " = " << K;
        err.info(ker.str());
      }
  }


namespace get_normalization_factor_impl
  {

    GiNaC::ex compute_normalization(const GiNaC::ex& term, service_locator& sl)
      {
        if(GiNaC::is_a<GiNaC::numeric>(term)) return term;

        if(GiNaC::is_a<GiNaC::symbol>(term))
          {
            // ignore if this symbol is the redshift z or is a parameter
            const auto& sym = GiNaC::ex_to<GiNaC::symbol>(term);

            const auto& sf = sl.get_symbol_factory();
            const auto& z = sf.get_z();

            if(static_cast<bool>(z == sym)) return GiNaC::numeric{1};
            if(sf.is_parameter(sym))
              {
                error_handler err;
                std::ostringstream msg;
                msg << WARNING_PARAMETER_APPEARS_IN_TIME_FUNCTION_RENORMALIZATION << ": '" << sym << "'";
                err.warn(msg.str());
              }

            return sym;
          }

        if(GiNaC::is_a<GiNaC::power>(term))
          {
            auto sub_norm = get_normalization_factor(term.op(0), sl);
            return GiNaC::pow(sub_norm, term.op(1));
          }

        return GiNaC::numeric{1};
      }

  }   // namespace get_normalization_factor_impl


GiNaC::ex get_normalization_factor(const time_function& tm, service_locator& sl)
  {
    // Want to canonicalize the time function.
    // This can be done by expanding it and collecting the numerical factors associated with the first term.
    // GiNaC guarantees that all expressions will come in a canonical order (even if we don't know what that
    // order is), so this will give consistent results for a single run of the code.
    // It won't necessarily give consistent answers *between* runs

    auto expr = tm.expand();

    using get_normalization_factor_impl::compute_normalization;

    // if expression is an add, extract numerical factors from its first term
    if(GiNaC::is_a<GiNaC::add>(expr)) return get_normalization_factor(expr.op(0), sl);

    // if expression is a mul, extract numerical factors
    if(GiNaC::is_a<GiNaC::mul>(expr))
      {
        GiNaC::ex norm{1};
        size_t ops = expr.nops();

        for(size_t i = 0; i < ops; ++i)
          {
            norm *= compute_normalization(expr.op(i), sl);
          }

        return norm;
      }

    // otherwise, determine normalization directly
    return compute_normalization(expr, sl);
  }


symmetrization_db build_symmetrizations(const GiNaC::ex& K, const initial_value_set& s)
  {
    // bin initial values by their symbol name
    // this gives us the individual groups over which we need to symmetrize
    using bin_db = std::map< GiNaC::symbol, std::vector<GiNaC::symbol> >;
    bin_db s_groups;

    for(auto t = s.value_cbegin(); t != s.value_cend(); ++t)
      {
        s_groups[t->get_symbol()].emplace_back(t->get_momentum());
      }

    // build atomic symmetrization database consisting of the identity transformation, ie. an empty exchange list
    symmetrization_db db = { GiNaC::exmap{} };

    // work through bins, adding elements to the symmetrization database for any with multiple occupancy
    for(auto& group : s_groups)
      {
        // get list of momenta associated with this symbol
        auto& original = group.second;
        auto group_size = original.size();

        // no need to symmetrize if just one instance
        if(group_size > 1)
          {
            // sort momentum list and copy to a separate list for generating permutations; need an explicit comparator
            // since operator< applied to
            std::sort(original.begin(), original.end(), std::less<GiNaC::symbol>{});
            auto perms = original;

            symmetrization_db new_db;
            do
              {
                // add new symmetrizations to each existing element in db and cache in new_db
                // later, we replace new_db with db
                for(const auto& perm : db)
                  {
                    GiNaC::exmap subs_map = perm;

                    for(unsigned int i = 0; i < perms.size(); ++i)
                      {
                        if(perms[i] != original[i]) subs_map[original[i]] = perms[i];
                      }

                    new_db.emplace_back(subs_map);
                  }
              }
            while(std::next_permutation(perms.begin(), perms.end(), std::less<GiNaC::symbol>{}));

            // replace db with new_db
            db.swap(new_db);
          }
      }

    return db;
  }


enum class factor_type { time, K, either };

factor_type classify_factor(const GiNaC::ex& expr, service_locator& loc)
  {
    // numbers can go anywhere
    if(GiNaC::is_a<GiNaC::numeric>(expr)) return factor_type::either;

    if(GiNaC::is_a<GiNaC::symbol>(expr))
      {
        const auto& sym = GiNaC::ex_to<GiNaC::symbol>(expr);

        auto& sf = loc.get_symbol_factory();
        const auto& z = sf.get_z();

        if(static_cast<bool>(sym == z)) return factor_type::time;
        if(sf.is_parameter(sym)) return factor_type::K;
        return factor_type::K;    // should perhaps be either, but other symbols are likely to be momenta and not allowed in time factor
      }

    if(GiNaC::is_a<GiNaC::add>(expr) || GiNaC::is_a<GiNaC::mul>(expr) || GiNaC::is_a<GiNaC::function>(expr) || GiNaC::is_a<GiNaC::power>(expr))
      {
        unsigned int t_only = 0;
        unsigned int K_only = 0;

        for(size_t i = 0; i < expr.nops(); ++i)
          {
            const auto& term = expr.op(i);
            auto type = classify_factor(term, loc);

            if(type == factor_type::time) ++t_only;
            if(type == factor_type::K) ++K_only;
          }

        if(t_only > 0 && K_only > 0)
          {
            std::ostringstream msg;
            msg << ERROR_INCONSISTENT_KERNEL_FACTOR << " " << expr;
            throw exception(msg.str(), exception_code::kernel_error);
          }

        if(t_only > 0) return factor_type::time;
        if(K_only > 0) return factor_type::K;
        return factor_type::either;
      }

    return factor_type::either;
  }


std::pair<GiNaC::ex, GiNaC::ex> sort_factors(const GiNaC::mul& expr, service_locator& loc)
  {
    GiNaC::ex tf{1};
    GiNaC::ex Kf{1};

    for(size_t i = 0; i < expr.nops(); ++i)
      {
        const auto& factor = expr.op(i);

        auto type = classify_factor(factor, loc);

        switch(type)
          {
            case factor_type::time:
              {
                tf *= factor;
                break;
              }

            case factor_type::K:
              {
                Kf *= factor;
                break;
              }

            case factor_type::either:
              {
                // prefer to move to kernel
                Kf *= factor;
                break;
              }
          }
      }

    return std::make_pair(tf, Kf);
  };


std::pair<GiNaC::ex, GiNaC::ex> partition_factor(const GiNaC::ex& expr, service_locator& loc)
  {
    GiNaC::ex tf{1};
    GiNaC::ex Kf{1};

    auto work_expr = GiNaC::collect_common_factors(expr.expand());

    if(GiNaC::is_a<GiNaC::mul>(work_expr))    // most common case, expr is a product of factors which we must sort
      {
        return sort_factors(GiNaC::ex_to<GiNaC::mul>(work_expr), loc);
      }

    // otherwise, factor is of a single type (even if it is a sum)
    auto type = classify_factor(work_expr, loc);

    switch(type)
      {
        case factor_type::time:
          {
            tf *= work_expr;
            break;
          }

        case factor_type::K:
          {
            Kf *= work_expr;
            break;
          }

        case factor_type::either:
          {
            // prefer to move to kernel
            Kf *= work_expr;
            break;
          }
      }

    return std::make_pair(tf, Kf);
  }
