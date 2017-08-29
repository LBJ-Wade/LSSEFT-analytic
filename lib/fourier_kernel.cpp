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
#include "utilities/GiNaC_utils.h"


namespace fourier_kernel_impl
  {
    
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
    
    
    void kernel::compose_substitution_rules(const kernel& rhs, GiNaC::exmap& mma_map)
      {
        // to do this, first check whether any of their substitution rules already exist in ours
        // if they do, add a suitable relabelling rule to mma_map
        // otherwise, insert a new rule provided there is no symbol collision
        for(const auto& v : rhs.vs)
          {
            const GiNaC::ex& label = v.first;
            const GiNaC::ex& value = v.second.subs(mma_map);
        
            const auto& label_symbol = GiNaC::ex_to<GiNaC::symbol>(label);
        
            // search for an existing remap rule which matches this one
            auto u = std::find_if(this->vs.begin(), this->vs.end(),
                                  [&](const subs_list::value_type& a) -> bool
                                    { return static_cast<bool>(a.second == value); });
        
            // if one exists, just add a relabelling rule
            if(u != this->vs.end())
              {
                mma_map[label_symbol] = u->first;
                break;
              }
        
            // otherwise need to add a new relabelling rule
        
            // first, is there a symbol collision?
            bool collision = false;
            const auto our_syms = this->iv.get_momenta();
            if(our_syms.find(label_symbol) != our_syms.end()) collision = true;
            if(this->vs.find(label_symbol) != this->vs.end()) collision = true;
        
            // if no collision, just keep the same symbol to avoid proliferation
            if(!collision)
              {
                this->vs[label_symbol] = value;
                break;
              }
        
            // otherwise, need to manufacture a new symbol
            auto relabel = this->sf.make_unique_momentum();
        
            mma_map[label_symbol] = relabel;
            this->vs[relabel] = value;
          }
      }
    
    
    kernel& kernel::operator+=(const kernel& rhs)
      {
        // we need to relabel momenta in the right-hand-side if they do not match
        // our current labelling
        auto our_mma = this->get_ordered_momenta();
        auto their_mma = rhs.get_ordered_momenta();
        
        // check that momenta are compatible in the sense that their (ordered) symbol names agree
        if(!std::equal(our_mma.cbegin(), our_mma.cend(),
                       their_mma.cbegin(), their_mma.cend(),
                       [](const momenta_list::value_type& a, const momenta_list::value_type& b) -> bool
                         { return a.get().get_symbol().get_name() == b.get().get_symbol().get_name(); }))
          throw exception(ERROR_KERNEL_INITIAL_VALUES_DISAGREE, exception_code::kernel_error);
        
        // build a substitution map for those momenta that disagree
        GiNaC::exmap mma_map;
        
        for(auto u = our_mma.cbegin(), v = their_mma.cbegin(); u != our_mma.cend() && v != their_mma.cend(); ++u, ++v)
          {
            const GiNaC::symbol& our_q = u->get().get_momentum();
            const GiNaC::symbol& their_q = u->get().get_momentum();

            // if the momentum labels disagree, then remap theirs
            if(our_q != their_q) mma_map[their_q] = our_q;
          }
        
        // also need to remap and merge any substitution rules in the right-hand side
        this->compose_substitution_rules(rhs, mma_map);
        
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
    
        // merge RHS substitution list with ours
        this->compose_substitution_rules(rhs, mma_map);
    
        // merge RHS initial value list with ours
        for(auto t = rhs.iv.value_cbegin(); t != rhs.iv.value_cend(); ++t)
          {
            const auto& q = t->get_momentum();
            auto u = mma_map.find(q);
            if(u != mma_map.end())
              {
                const auto& u_sym = GiNaC::ex_to<GiNaC::symbol>(u->second);
                this->iv.insert(t->relabel_momentum(u_sym));
              }
            else
              {
                this->iv.insert(*t);
              }
          }

        // build final expression, performing any necessary index (or other) relabellings on RHS
        using detail::relabel_index_product;
        auto temp = simplify_index(relabel_index_product(this->K, rhs.K.subs(mma_map), this->sf));
        this->K = temp;

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
    
    
    kernel operator*(const GiNaC::ex a, const kernel& b)
      {
        kernel c = b;
        c.tm *= a;
        return c;
      }
    
    
    kernel operator*(const kernel& a, const kernel& b)
      {
        kernel c = a;
        c *= b;
        return c;
      }
    
    
    kernel diff_z(const kernel& a)
      {
        kernel b = a;
        const auto& z = a.sf.get_z();
        b.tm = GiNaC::diff(b.tm, z);
        
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
        this->K *= f;
        return *this;
      }
    
    
    kernel& kernel::multiply_kernel(GiNaC::ex f, GiNaC::symbol s, GiNaC::ex rule)
      {
        const auto& our_syms = this->iv.get_momenta();
        const auto& eps = this->sf.get_regulator();

        // ensure that f only involves symbols in the initial value set or s (or the regulator epsilon)
        auto expr_syms = get_expr_symbols(f);
    
        // remove regulator epsilon from the used set
        auto t = expr_syms.find(eps);
        if(t != expr_syms.end()) expr_syms.erase(t);

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
        
        this->K *= f;
        this->vs[s] = rule;
        
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
  

void validate_subslist(const initial_value_set& s, const fourier_kernel_impl::subs_list& vs)
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


void validate_momenta(const initial_value_set& s, const fourier_kernel_impl::subs_list& vs, const GiNaC::ex& K,
                      const GiNaC::symbol& eps, bool silent)
  {
    using fourier_kernel_impl::subs_list;
    
    auto used = get_expr_symbols(K);
    auto avail = s.get_momenta();
    
    // insert any symbols from the substitution list vs into the available set
    std::for_each(vs.begin(), vs.end(), [&](const subs_list::value_type& v) -> void
      {
        auto sym = GiNaC::ex_to<GiNaC::symbol>(v.first);
        avail.insert(sym);
      });
    
    GiNaC_symbol_set used_not_avail;
    GiNaC_symbol_set avail_not_used;
    
    // remove regulator epsilon from the used set
    auto t = used.find(eps);
    if(t != used.end()) used.erase(t);
    
    // perform set differencing to find mismatch between available and used symbols
    std::set_difference(used.cbegin(), used.cend(), avail.cbegin(), avail.cend(),
                        std::inserter(used_not_avail, used_not_avail.end()), std::less<GiNaC::symbol>{});
    std::set_difference(avail.cbegin(), avail.cend(), used.cbegin(), used.cend(),
                        std::inserter(avail_not_used, avail_not_used.end()), std::less<GiNaC::symbol>{});
    
    // raise exception on attempt to use a Fourier momentum that isn't available
    if(!used_not_avail.empty())
      {
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
