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


    class Wick_data
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! constructor captures data
        Wick_data(GiNaC::ex Pks_);
        
        //! destructor is default
        ~Wick_data() = default;


        // INTERFACE

      public:

        //! extract Wick product string
        const GiNaC::ex& get_Wick_string() const { return this->Pk_string; }
        
        
        // INTERNAL DATA
        
      private:
        
        //! substitution map for kernel 1
        GiNaC::exmap ker1_subs;
        
        //! substitution map for kernel 2
        GiNaC::exmap ker2_subs;
        
        //! list of remaining free (ie. loop) momenta
        GiNaC_symbol_set loop_momenta;
        
        //! GiNaC expression representing the string of correlation functions
        //! produced by this set of Wick contraction
        GiNaC::ex Pk_string;

      };
    
    
    //! the set of all possible contractions is a list of contraction_data objects
    using Wick_set = std::vector< std::unique_ptr<Wick_data> >;

    
    //! contractions builds a list of contractions between initial value strings
    class contractions
      {
        
        // TYPES
        
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
        
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! constructor takes a list of initial_value sets,
        //! a list of corresponding external momenta,
        //! and constructs a set of all possible contractions between them
        template <size_t N>
        contractions(std::array<initial_value_set, N> groups, std::array<GiNaC::ex, N> kext);
        
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

        //! build the data needed to construct a Wick product from a contraction group
        template <size_t N>
        void build_Wick_product(const contraction_group& gp,
                                const std::array<initial_value_set, N>& groups,
                                const std::array<GiNaC::ex, N>& kext);

        
        // INTERNAL DATA
        
      private:

        //! set of Wick contractions
        std::unique_ptr<Wick_set> items;
      
      };


    template <size_t N>
    contractions::contractions(std::array<initial_value_set, N> groups, std::array<GiNaC::ex, N> kext)
      : items(std::make_unique<Wick_set>())
      {
        // total number of fields should be even since we currently include only Gaussian contractions
        // that pair together exactly two fields
        size_t num = 0;
        for(const auto& group : groups)
          {
            num += group.size();
          }

        if(num % 2 != 0)
          throw exception(ERROR_ODD_CONTRACTIONS, exception_code::contraction_error);

        // build a single list of all fields; we'll use this to pull out contractions pair-by-pair
        iv_list ivs;
        ivs.reserve(num);

        for(size_t i = 0; i < N; ++i)
          {
            const auto& group = groups[i];
            for(auto t = group.value_cbegin(); t != group.value_cend(); ++t)
              {
                // emplace an iterator to this initial field, tagged with the id of the group
                // from which it comes
                ivs.emplace_back(t, i);
              }
          }

        // now build all possible Wick pairs, keeping only those that produce connected
        // correlation functions
        auto ctrs = this->enumerate_contractions(num, ivs);

        // iterate over contractions groups held in ctrs
        for(const auto& gpp : *ctrs)
          {
            const auto& gp = *gpp;

            // convert this contraction group into the data we need to supply -- eg. GiNaC substitution lists,
            // lists of external momenta, etc ...
            this->build_Wick_product(gp, groups, kext);
          }
      }


    template <size_t N>
    void contractions::build_Wick_product(const contraction_group& gp,
                                          const std::array<initial_value_set, N>& groups,
                                          const std::array<GiNaC::ex, N>& kext)
      {
        // need to build a string of power spectra representing the Wick product in gp
        Pk_string Ps;

        // keep track of which momenta are integrated over, corresponding to loops
        GiNaC_symbol_set loop_momenta;

        for(const auto& prod : gp)
          {
            const auto& left = prod.first;
            const auto& right = prod.second;

            const auto& left_sym = left.first->get_symbol();
            const auto& left_mom = left.first->get_momentum();

            const auto& right_sym = right.first->get_symbol();
            const auto& right_mom = right.first->get_momentum();

            // place symbols into canonical order, inherited from std::less<> applied to GiNaC
            // symbols (recall we define this ourselves to give lexical order on the symbol names)
            GiNaC::symbol l;
            GiNaC::symbol r;
            if(std::less<>{}(left_sym, right_sym))
              {
                l = left_sym;
                r = right_sym;
              }
            else
              {
                l = right_sym;
                r = left_sym;
              }

            // take momentum from first field (the choice is arbitrary)
            Ps.emplace_back(cfs::Pk(l, r, left_mom), left.second);
          }

        // convert Ps to a GiNaC product
        GiNaC::ex W = 1;
        for(const auto& factor : Ps)
          {
            W *= factor.first;
          }

        this->items->emplace_back(std::make_unique<Wick_data>(W));
      }


  }   // namespace detail


#endif //LSSEFT_ANALYTIC_CONTRACTIONS_H
