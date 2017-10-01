//
// Created by David Seery on 25/08/2017.
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

#ifndef LSSEFT_ANALYTIC_CONTRACTIONS_H
#define LSSEFT_ANALYTIC_CONTRACTIONS_H


#include <array>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <memory>

#include "lib/initial_value.h"

#include "utilities/hash_combine.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"

#include "relabel_product.h"

#include "ginac/ginac.h"


namespace cfs
  {

    //! declare a GiNaC function to represent a power spectrum
    //! this is really just a placeholder function whose job is to accept three arguments:
    //! the first two specify the correlation function in question and the third specifies
    //! the momentum
    DECLARE_FUNCTION_3P(Pk)

  }   // namespace cfs = correlation functions



namespace detail
  {

    //! graph represents the Feynman graph built by performing Wick contractions;
    //! we need it in order to know whether we have produced a connected or disconnected
    //! correlation function
    class graph
      {

        // TYPES

      protected:

        //! represent a graph as an list of edges source -> destination,
        //! where source and destination label vertices
        //! this is used to determine whether a given Wick contraction produces
        //! a connected or disconnected graph
        //! we store the edges as a map from source vertex to the list of destination
        //! vertices available from that source
        using edge_db = std::map< size_t, std::set<size_t> >;

        //! represent list of vertices
        using vertex_db = std::set<size_t>;


        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor generates an empty graph
        graph() = default;

        //! destructor is default
        ~graph() = default;


        // EDGE HANDLING

      public:

        //! add an edge from source -> destination
        //! (we track orientation of edges but this information isn't really used)
        void add_edge(size_t source, size_t dest);


        // QUERY

      public:

        //! determine whether the graph is connected
        bool is_connected() const;


        // INTERNAL DATA

      private:

        //! storage for the graph
        edge_db edges;

        //! storage for vertices
        vertex_db vertices;

      };
    
    
    //! type for substitution rules for Rayleigh momenta
    using subs_map = GiNaC::exmap;


    class Wick_data
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! constructor captures data
        template <typename RemapContainer>
        Wick_data(GiNaC::ex Pks_, RemapContainer& mma_, GiNaC_symbol_set lm_);
        
        //! destructor is default
        ~Wick_data() = default;


        // INTERFACE

      public:

        //! extract Wick product string
        const GiNaC::ex& get_Wick_string() const { return this->Pk_string; }
        
        //! extract set of loop momenta
        const GiNaC_symbol_set& get_loop_momenta() const { return this->loop_momenta; }
        
        //! extract substitution rules for each kernel
        const std::vector<subs_map>& get_substitution_rules() const { return this->mma_map; }
        
        
        // INTERNAL DATA
        
      private:
        
        //! substitution rules for each kernel
        std::vector<subs_map> mma_map;
        
        //! list of remaining free (ie. loop) momenta
        GiNaC_symbol_set loop_momenta;
        
        //! GiNaC expression representing the string of correlation functions
        //! produced by this set of Wick contraction
        GiNaC::ex Pk_string;

      };
    
    
    template <typename RemapContainer>
    Wick_data::Wick_data(GiNaC::ex Pks_, RemapContainer& mma_, GiNaC_symbol_set lm_)
      : Pk_string(std::move(Pks_)),
        loop_momenta(std::move(lm_))
      {
        // copy contents of mma into remap rules
        std::move(mma_.begin(), mma_.end(), std::back_inserter(this->mma_map));
      }
    
    
    //! the set of all possible contractions is a list of contraction_data objects
    using Wick_set = std::vector< std::unique_ptr<Wick_data> >;

    
    //! contractions builds a list of contractions between initial value strings
    class contractions
      {
        
        // TYPES
        
      public:
        
        //! initial value groups
        template <size_t N>
        using iv_group = std::array<initial_value_set, N>;
        
        //! group of external momenta
        template <size_t N>
        using kext_group = std::array<GiNaC::ex, N>;
        
      protected:

        //! iterator to an initial_value item
        using iv_item = initial_value_set::iv_db::const_iterator;

        //! a contraction element specifies a single initial value item and a label to indicate which group (=vertex) it belongs to
        using iv_element = std::pair<iv_item, size_t>;

        //! list of initial value items
        using iv_list = std::vector<iv_element>;

        //! a contraction is a pairing between initial_value items
        using contraction = std::pair<iv_element, iv_element>;
        
        //! a contraction group of a list of contractions
        using contraction_group = std::vector<contraction>;
        
        //! a contraction set of a list of contraction groups (usually representing all possible contractions of an index set)
        using contraction_set = std::vector< std::unique_ptr<contraction_group> >;

        //! an element of a string of power spectra generated by Wick products
        using Pk_element = std::pair<GiNaC::ex, size_t>;

        //! a string of such power spectrum elements
        using Pk_string = std::vector<Pk_element>;
    
    
        //! array of exmaps used to remap kernel momentum labels after assigning loop momenta
        template <size_t N>
        using remap_group = std::array<GiNaC::exmap, N>;
    
        //! convenience type holding a std::set of momentum labels that are yet to be fixed (during assignment of loop momenta)
        template <size_t N>
        using symbol_group = std::array<GiNaC_symbol_set, N>;
    
    
    
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! constructor takes a list of initial_value sets ('clusters'),
        //! a list of corresponding external momenta,
        //! and constructs a set of all possible contractions between them
        template <size_t N>
        contractions(iv_group<N> clusters, kext_group<N> kext, service_locator& lc_);
        
        //! destructor is default
        ~contractions() = default;


        // SERVICES

      public:

        //! get data for available Wick contractions
        const Wick_set& get() const { return *this->items; }


        // INTERNAL API

      protected:

        //! build and enumerate all possible contractions for a pair of index sets
        std::unique_ptr<contraction_set>
        enumerate_contractions(size_t num, const iv_list& ivs) const;
        
        //! compute the loop order of a given contraction group
        template <size_t N>
        size_t loop_order(const contraction_group& gp, const iv_group<N>& clusters);
        
        //! assign loop momenta
        template <size_t N, typename LabelMap, typename UnassignedGroup>
        void assign_loop_momenta(const contraction_group& gp, const iv_group <N>& clusters, LabelMap& mma_map,
                                 UnassignedGroup& unassigned, GiNaC_symbol_set& loop_momenta);
        
        //! assign external momenta
        template <size_t N, typename LabelMap, typename UnassignedGroup>
        void assign_external_momenta(const iv_group<N>& clusters, const kext_group<N>& kext, LabelMap& mma_map,
                                     UnassignedGroup& unassigned);
        
        //! evaluate all Wick contractions
        template <size_t N, typename LabelMap, typename UnassignedGroup>
        void evaluate_Wick_contractions(const contraction_group& gp, const iv_group<N>& clusters,
                                        LabelMap& mma_map, const UnassignedGroup& unassigned, Pk_string& Ps);
        
        //! build the data needed to construct a Wick product from a contraction group
        template <size_t N>
        void build_Wick_product(const contraction_group& gp, const iv_group<N>& clusters, const kext_group<N>& kext);

        
        // INTERNAL DATA
        
      private:
        
        //! cache reference to service locator
        service_locator& loc;

        //! set of Wick contractions
        std::unique_ptr<Wick_set> items;
      
      };


    template <size_t N>
    contractions::contractions(iv_group<N> clusters, kext_group<N> kext, service_locator& lc_)
      : items(std::make_unique<Wick_set>()),
        loc(lc_)
      {
        // total number of fields should be even since we currently include only Gaussian contractions
        // that pair together exactly two fields
        size_t num = 0;
        for(const auto& cluster : clusters)
          {
            num += cluster.size();
          }

        if(num % 2 != 0)
          throw exception(ERROR_ODD_CONTRACTIONS, exception_code::contraction_error);

        // build a single list of all fields; we'll use this to pull out contractions pair-by-pair
        iv_list ivs;
        ivs.reserve(num);

        for(size_t i = 0; i < N; ++i)
          {
            const auto& cluster = clusters[i];
            for(auto t = cluster.value_cbegin(); t != cluster.value_cend(); ++t)
              {
                // emplace an iterator to this initial field, tagged with the id of the cluster
                // from which it comes
                ivs.emplace_back(t, i);
              }
          }

        // now build all possible Wick pairs, keeping only those that produce connected
        // correlation functions; we check this by building a representation of the
        // Wick graph and checking it for connectedness
        auto ctrs = this->enumerate_contractions(num, ivs);

        // iterate over contractions groups held in ctrs
        for(const auto& gpp : *ctrs)
          {
            const auto& gp = *gpp;

            // convert this contraction group into the data we need to supply -- eg. GiNaC substitution lists,
            // lists of external momenta, etc ...
            this->build_Wick_product(gp, clusters, kext);
          }
      }
    
    
    template <size_t N, typename LabelMap, typename UnassignedGroup>
    void
    contractions::assign_external_momenta(const iv_group<N>& clusters, const kext_group<N>& kext, LabelMap& mma_map,
                                          UnassignedGroup& unassigned)
      {
        for(size_t i = 0; i < N; ++i)
          {
            const auto& labels = unassigned[i];
            
            if(labels.size() < 1)
              throw exception(ERROR_COULD_NOT_ASSIGN_EXTERNAL_MOMENTUM, exception_code::contraction_error);
            
            // if there is *more than one* unassigned momentum left in this cluster, then do nothing
            // currently, we assume all these free fields will be attached to an external momentum in another cluster
            // TODO: this may need revisiting in future, but currently seems OK at least for the 1-loop P and B
            if(labels.size() > 2) continue;
            
            // we are guaranteed there is exactly a single unassigned momentum left, so pick it out
            const auto label = *labels.begin();
            
            // Each cluster has a momentum conservation delta-function of the form delta(k_ext - q1 - q2 - ... - qn).
            // At this stage, our assumptions guarantee that all but one of the qi are now fixed
            // The remaining one, say qj, can be written qj = k_exit - q1 - ... - qn
            const auto& ks = clusters[i].get_momenta();
            GiNaC::ex sum = kext[i];
            for(const auto& k : ks)
              {
                if(k != label) sum -= k;
              }

            // the q_i have all been assigned, so we should make the appropriate replacements
            // to get their final values
            auto relabelled_sum = sum.subs(mma_map[i]);
            
            // store this final value in the replacement map
            mma_map[i][label] = relabelled_sum;
            
            // remove this label from the unassigned list
            unassigned[i].erase(unassigned[i].begin());
          }
      }
    
    
    template <size_t N, typename LabelMap, typename UnassignedGroup>
    void
    contractions::assign_loop_momenta(const contraction_group& gp, const iv_group<N>& clusters, LabelMap& mma_map,
                                      UnassignedGroup& unassigned, GiNaC_symbol_set& loop_momenta)
      {
        auto& sf = this->loc.get_symbol_factory();

        size_t loop_order = this->loop_order(gp, clusters);
    
        // take a working copy of the list of Wick contractions
        std::list< contraction > gp_copy;
        std::copy(gp.begin(), gp.end(), std::back_inserter(gp_copy));
    
        // predicate to search for candidate Wick contractions
        auto find_candidate = [&](const contraction& c) -> bool
          {
            // work out which clusters the legs of this contraction belong to
            const auto& clust1 = c.first.second;
            const auto& clust2 = c.second.second;
        
            // if there is only a single momentum remaining unassigned at either vertex, then this
            // contraction isn't a loop -- it must communicate with an external momentum
            if(unassigned[clust1].size() <= 1) return false;
            if(unassigned[clust2].size() <= 1) return false;
        
            // otherwise, we can assign a loop momentum
            return true;
          };
    
        // attempt to assign loop momentum labels
        for(size_t i = 0; i < loop_order; ++i)
          {
            // search gp_copy for a candidate contraction to assign as a loop
            auto t = std::find_if(gp_copy.begin(), gp_copy.end(), find_candidate);
        
            if(t == gp_copy.end())
              throw exception(ERROR_COULD_NOT_ASSIGN_LOOP_MOMENTUM, exception_code::contraction_error);
        
            // remove these momentum labels from the unassigned list
            const auto& clust1 = t->first.second;
            const auto& clust2 = t->second.second;
        
            const auto& iv1 = *t->first.first;
            const auto& iv2 = *t->second.first;
        
            auto t1 = unassigned[clust1].find(iv1.get_momentum());
            auto t2 = unassigned[clust2].find(iv2.get_momentum());
            
            if(t1 == unassigned[clust1].end() || t2 == unassigned[clust2].end())
              throw exception(ERROR_COULD_NOT_ASSIGN_LOOP_MOMENTUM, exception_code::contraction_error);
        
            unassigned[clust1].erase(t1);
            unassigned[clust2].erase(t2);
            
            // manufacture a new loop momentum variable
            auto l = sf.make_unique_loop_momentum();

            // generate replacement rules for the momenta we have paired up
            mma_map[clust1][iv1.get_momentum()] = l;
            mma_map[clust2][iv2.get_momentum()] = -l;
            
            // insert the label l in the list of loop momenta
            loop_momenta.insert(l);
        
            // remove this contraction from the candidate list
            gp_copy.erase(t);
          }
      }
    
    
    template <size_t N, typename LabelMap, typename UnassignedGroup>
    void
    contractions::evaluate_Wick_contractions(const contraction_group& gp, const iv_group<N>& clusters,
                                             LabelMap& mma_map, const UnassignedGroup& unassigned, Pk_string& Ps)
      {
        for(const auto& prod : gp)
          {
            const auto& left = prod.first;
            const auto& right = prod.second;
            
            const auto& left_clust = left.second;
            const auto& right_clust = right.second;
        
            const auto& left_sym = left.first->get_symbol();
            const auto& left_mom = left.first->get_momentum();
        
            const auto& right_sym = right.first->get_symbol();
            const auto& right_mom = right.first->get_momentum();
    
            // first, if both momenta appear in the 'unassigned' set then this implementation doesn't
            // currently know what to do. Give up in despair.
            bool left_unassigned = unassigned[left_clust].find(left_sym) != unassigned[left_clust].end();
            bool right_unassigned = unassigned[right_clust].find(right_sym) != unassigned[right_clust].end();
            
            if(left_unassigned && right_unassigned)
              throw exception(ERROR_COULD_NOT_EVALUATE_WICK_CONTRACTION, exception_code::contraction_error);
    
            // place symbols into canonical order, inherited from std::less<> applied to GiNaC
            // symbols (recall we define this ourselves to give lexical order on the symbol names)
            GiNaC::symbol l = left_sym;
            GiNaC::symbol r = right_sym;
            if(std::less<>{}(right_sym, left_sym)) std::swap(l,r);
            
            // we have a choice which momentum to use
            // it's preferable to use either a simple external momentum or a simple loop momentum
            GiNaC::ex q;
            
            // if either the left- or right-hand momentum is unassigned, we always use the other one
            if(left_unassigned)
              {
                Ps.emplace_back(cfs::Pk(l, r, right_mom).subs(mma_map[right_clust]), right_clust);
                
                // add a rewriting rule for this contraction
                mma_map[left_clust][left_sym] = -right_mom.subs(mma_map[right_clust]);
                break;
              }
            if(right_unassigned)
              {
                Ps.emplace_back(cfs::Pk(l, r, left_mom).subs(mma_map[left_clust]), left_clust);
                
                // add a rewriting rule for this contraction
                mma_map[right_clust][right_sym] = -left_mom.subs(mma_map[left_clust]);
                break;
              }
            
            // otherwise, if the left-momentum is a simple symbol then we should use it;
            // this will catch cases where the momentum is exactly an external momentum or a simple loop
            auto t1 = mma_map[left_clust].find(left_sym);
            if(t1 != mma_map[left_clust].end() && GiNaC::is_a<GiNaC::symbol>(t1->second))
              {
                Ps.emplace_back(cfs::Pk(l, r, left_mom).subs(mma_map[left_clust]), left_clust);
                break;
              }
            
            // same for right-momentum
            auto t2 = mma_map[right_clust].find(right_sym);
            if(t2 != mma_map[right_clust].end() && GiNaC::is_a<GiNaC::symbol>(t2->second))
              {
                Ps.emplace_back(cfs::Pk(l, r, right_mom).subs(mma_map[right_clust]), right_clust);
                break;
              }
        
            // nothing to choose between the LHS and RHS momenta, so just pick one
            Ps.emplace_back(cfs::Pk(l, r, left_mom).subs(mma_map[left_clust]), left_clust);
          }
      }


    template <size_t N>
    void
    contractions::build_Wick_product(const contraction_group& gp, const iv_group<N>& clusters, const kext_group<N>& kext)
      {
        // need to build a string of power spectra representing the Wick product in gp
        Pk_string Ps;

        // keep track of which momenta are integrated over, corresponding to loops
        GiNaC_symbol_set loop_momenta;
        
        // keep track of which momenta need Rayleigh expansion
        subs_map Rayleigh_momenta;
        
        // STEP 1. Set up replacement rules for relabelling the momenta carried by each kernel
        remap_group<N> mma_map;
        symbol_group<N> unassigned;

        // populate 'unassigned' with the momentum labels for each group
        for(size_t i = 0; i < N; ++i)
          {
            unassigned[i] = clusters[i].get_momenta();
          }
        
        this->assign_loop_momenta(gp, clusters, mma_map, unassigned, loop_momenta);
        this->assign_external_momenta(clusters, kext, mma_map, unassigned);
        
        // STEP 2. Construct Wick contractions
        this->evaluate_Wick_contractions(gp, clusters, mma_map, unassigned, Ps);
        
        // STEP 3. Convert Ps to a GiNaC product corresponding to the outcome of all Wick contractions
        GiNaC::ex W = 1;
        std::for_each(Ps.begin(), Ps.end(), [&](const auto& f) -> void { W *= f.first; });

        // STEP 4. Store all this data
        this->items->emplace_back(std::make_unique<Wick_data>(W, mma_map, loop_momenta));
      }
    
    
    template <size_t N>
    size_t
    contractions::loop_order(const contraction_group& gp, const iv_group<N>& clusters)
      {
        // determine number of contractions
        size_t contractions = gp.size();
        
        // determine number of initial field in the full cluster size
        size_t fields = 0;
        for(const auto& cluster : clusters)
          {
            fields += cluster.size();
          }
    
        // ensure #fields is a multiple of 2
        if(fields % 2 != 0)
          throw exception(ERROR_ODD_CONTRACTIONS, exception_code::contraction_error);
    
        // there are #fields independent momentum integrals
        // each contraction generates one momentum-conservation delta function
        // there are also #clusters momentum-conservation delta functions, one for each cluster
        // finally, one of these delta functions factorizes out to give global conservation of momentum
        
        // so, the number of unconstrained momentum integrations is
        // #fields - #clusters - #contractions + 1
        return fields - clusters.size() - contractions + 1;
      }
    
    
  }   // namespace detail


#endif //LSSEFT_ANALYTIC_CONTRACTIONS_H
