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

#include "fourier_kernel.h"
#include "detail/relabel_product.h"
#include "detail/Rayleigh_momenta.h"
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
        std::hash<std::string> string_hasher;

        // to hash the time expression, expand it completely and print
        std::ostringstream time_string;
        time_string << this->tm.expand();

        size_t h = string_hasher(time_string.str());

        // to hash the initial value set, order its symbols lexicographically
        const auto symbols = this->get_ordered_iv_symbols();

        std::string symbol_string;
        std::for_each(symbols.begin(), symbols.end(),
                      [&](const GiNaC::symbol& e) -> std::string
                        { return symbol_string += e.get_name(); });

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


    kernel::kernel(GiNaC::ex K_, initial_value_set iv_, time_function tm_, subs_list vs_, symbol_factory& sf_)
      : K(std::move(K_)),
        tm(std::move(tm_)),
        iv(std::move(iv_)),
        vs(std::move(vs_)),
        sf(sf_)
      {
        // normalize the time function, redistributing factors into the kernel if needed
        auto norm = get_normalization_factor(tm, sf);
        tm /= norm;
        K *= norm;
      }


    kernel::kernel(initial_value_set iv_, symbol_factory& sf_)
      : K(1),
        tm(1),
        iv(std::move(iv_)),
        sf(sf_)
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
        using detail::merge_Rayleigh_rules;
        auto relabel_map = merge_Rayleigh_rules(this->vs, rhs.vs, this->iv.get_momenta(), mma_map, this->sf);
        std::copy(relabel_map.begin(), relabel_map.end(), std::inserter(mma_map, mma_map.begin()));

        // now perform relabelling in the kernel
        // note there's no need to perform *index* relabelling since this is a sum -- relabelling indices
        // is needed only in a product
        auto temp = simplify_index(this->K + rhs.K.subs(mma_map));
        this->K = temp;

        return *this;
      }
    
    
    kernel& kernel::operator*=(const kernel& rhs)
      {
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
            auto relabel = this->sf.make_unique_momentum();
            
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
        using detail::merge_Rayleigh_rules;
        auto relabel_map = merge_Rayleigh_rules(this->vs, rhs.vs, this->iv.get_momenta(), mma_map, this->sf);
        std::copy(relabel_map.begin(), relabel_map.end(), std::inserter(mma_map, mma_map.begin()));

        // build final expression, performing any necessary index (or other) relabellings on RHS
        using detail::relabel_index_product;
        auto temp = simplify_index(relabel_index_product(this->K, rhs.K.subs(mma_map), this->sf));
        this->K = temp;

        // renormalize time function
        auto norm = get_normalization_factor(tm, sf);
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

        // explicitly distribute a over a top-level sum if one is present
        // this is because GiNaC is very reluctant to simplify cancelling factors in such cases,
        // meaning that we get uncancelled factors of 1+z in the result, even though they all cancel down to 1
        if(GiNaC::is_a<GiNaC::add>(c.tm))
          {
            GiNaC::ex temp{0};
            const size_t nops = c.tm.nops();

            for(size_t i = 0; i < nops; ++i)
              {
                temp += c.tm.op(i).expand() * a;
              }

            c.tm = GiNaC::collect_common_factors(temp);
          }
        else
          {
            auto temp = GiNaC::collect_common_factors(c.tm.expand() * a);
            c.tm = temp;
          }


        // renormalize time function
        auto norm = get_normalization_factor(c.tm, c.sf);
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
        kernel b = a;
        const auto& z = a.sf.get_z();
        b.tm = GiNaC::diff(b.tm, z);

        // renormalize time function
        auto norm = get_normalization_factor(b.tm, b.sf);
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
        const auto& params = this->sf.get_parameters();
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
        auto norm = get_normalization_factor(tm, sf);
        this->tm /= norm;
        this->K *= norm;

        return *this;
      }
    
    
    kernel& kernel::multiply_kernel(GiNaC::ex f, GiNaC::symbol s, GiNaC::ex rule)
      {
        auto our_syms = this->iv.get_momenta();
        const auto& params = this->sf.get_parameters();
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
        if(t == this->vs.end()) this->vs[s] = rule;

        // renormalize time function
        auto norm = get_normalization_factor(tm, sf);
        this->tm /= norm;
        this->K *= norm;

        return *this;
      }
    
    
    void kernel::write(std::ostream& out) const
      {
        out << "  time factor = " << this->tm << '\n';
    
        out << "  IC set =";
        for(auto u = this->iv.value_cbegin(); u != this->iv.value_cend(); ++u)
          {
            out << " " << u->get_symbol() << "(" << u->get_momentum() << ")";
          }
        out << '\n';
    
        if(!this->vs.empty())
          {
            out << "  substitution list = ";
            for(const auto& u : vs)
              {
                out << " " << u.first << " -> " << u.second << ";";
              }
            out << '\n';
          }
    
        out << "  kernel = " << this->K << '\n';
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
            msg << WARNING_KERNEL_EXPRESSION << " = " << K;
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
        ker << WARNING_KERNEL_EXPRESSION << " = " << K;
        err.info(ker.str());
      }
  }


namespace get_normalization_factor_impl
  {

    GiNaC::ex compute_normalization(const GiNaC::ex& term, symbol_factory& sf)
      {
        if(GiNaC::is_a<GiNaC::numeric>(term)) return term;

        if(GiNaC::is_a<GiNaC::symbol>(term))
          {
            // ignore if this symbol is the redshift z
            const auto& z = sf.get_z();
            if(static_cast<bool>(z == term)) return GiNaC::numeric{1};
            return term;
          }

        if(GiNaC::is_a<GiNaC::power>(term))
          {
            auto sub_norm = get_normalization_factor(term.op(0), sf);
            return GiNaC::pow(sub_norm, term.op(1));
          }

        return GiNaC::numeric{1};
      }

  }   // namespace get_normalization_factor_impl


GiNaC::ex get_normalization_factor(const time_function& tm, symbol_factory& sf)
  {
    // Want to canonicalize the time function.
    // This can be done by expanding it and collecting the numerical factors associated with the first term.
    // GiNaC guarantees that all expressions will come in a canonical order (even if we don't know what that
    // order is), so this will give consistent results for a single run of the code.
    // It won't necessarily give consistent answers *between* runs

    auto expr = tm.expand();

    using get_normalization_factor_impl::compute_normalization;

    // if expression is an add, extract numerical factors from its first term
    if(GiNaC::is_a<GiNaC::add>(expr)) return get_normalization_factor(expr.op(0), sf);

    // if expression is a mul, extract numerical factors
    if(GiNaC::is_a<GiNaC::mul>(expr))
      {
        GiNaC::ex norm{1};
        size_t ops = expr.nops();

        for(size_t i = 0; i < ops; ++i)
          {
            norm *= compute_normalization(expr.op(i), sf);
          }

        return norm;
      }

    // otherwise, determine normalization directly
    return compute_normalization(expr, sf);
  }
