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

#include "contractions.h"


namespace cfs
  {

    // register Pk as a dummy function
    REGISTER_FUNCTION(Pk, dummy())

  }   // namespace cfs = correlation functions


namespace detail
  {

    size_t double_factorial(size_t N)
      {
        size_t x = 1;
        for(size_t i = N; i > 1; i -= 2)
          {
            x *= i;
          }

        return x;
      }


    std::unique_ptr<contractions::contraction_set>
    contractions::enumerate_contractions(size_t num, const iv_list& ivs) const
      {
        // allocate a unique_ptr for the contraction set
        auto ctrs = std::make_unique<contraction_set>();

        // compute number of possible contractions = (N-1)!! where N is the size of the set
        size_t num_contractions = double_factorial(num-1);

        // work through total number of contractions, building up the pairings associated with each item
        for(size_t i = 0; i < num_contractions; ++i)
          {
            // allocate storage for this contraction group
            auto group = std::make_unique<contraction_group>();

            // copy list of initial values, which will be used to populate the contraction elements
            auto ivs_copy = ivs;

            // use an edge database to represent the Feynman graph we are effectively setting up;
            // we'll need to know whether it is connected
            detail::graph G;

            // divisor keeps track of numerical factor needed to strip out next term in the double factorial
            size_t divisor = num_contractions;
            size_t v = i;

            for(size_t j = num - 1; j >= 1; j -= 2)
              {
                divisor /= j;

                // compute which element we're going to pair with the first element remaining in ivs_copy;
                // note the +1, because we start counting from the *second* element in ivs_copy
                auto pair_elt = v / divisor + 1;

                // update for next pairing
                v = v % divisor;

                // extract iv_element objects for each end of the pairing
                const auto& p1 = ivs_copy[0];
                const auto& p2 = ivs_copy[pair_elt];

                // add an edge from p1 to p2
                G.add_edge(p1.second, p2.second);

                // store this pairing
                group->emplace_back(p1, p2);

                // remove this pair of fields from ivs_copy, ready for next contraction
                ivs_copy.erase(ivs_copy.begin() + pair_elt);
                ivs_copy.erase(ivs_copy.begin());

                // manually break out of for-loop if j=1; subtracting two will wrap around since
                // j is unsigned
                if(j == 1) break;
              }

            // if the graph we built was connected, store this contraction group
            if(G.is_connected()) ctrs->push_back(std::move(group));
          }

        return ctrs;
      }

    void contractions::build_Wick_product(const contractions::contraction_group& gp)
      {
        // need to build a string of power spectra representing the Wick product in gp
        Pk_string Ps;

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

            // take momentum from first field
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


    void graph::add_edge(size_t source, size_t dest)
      {
        // insert source and destination vertices if they have not already been encountered
        this->vertices.insert(source);
        this->vertices.insert(dest);

        // insert source -> dest edge
        auto& dest_set = this->edges[source];
        dest_set.insert(dest);
      }


    template <typename EdgeContainer, typename VertexContainer>
    void visit_all(const EdgeContainer& es, VertexContainer& vs, size_t v)
      {
        // if vs is already empty, return; we have already visited all nodes,
        // so there is no need to carry on
        if(vs.empty()) return;

        // remove vertex v from the container; this means we cannot recurse infinitely if the
        // graph contains a cycle
        auto t = vs.find(v);
        if(t == vs.end()) return;

        vs.erase(t);

        // find vertices that can be reached from this source; return if there are none -- this is a leaf node
        auto u = es.find(v);
        if(u == es.end()) return;

        for(const auto& d : u->second)
          {
            visit_all(es, vs, d);
          }
      }


    bool graph::is_connected() const
      {
        // define an empty graph to be disconnected
        if(this->edges.empty() || this->vertices.empty()) return false;

        // can assume there is at least one vertex and at least one line.
        // To determine whether the graph is connected we recursively visit all the vertices
        // attached to any single vertex in the vertex set.
        // If we visit all vertices then the graph is connected.
        // If not then it is disconnected
        vertex_db vs = this->vertices;

        visit_all(this->edges, vs, *vs.begin());

        return vs.empty();
      }


    Wick_data::Wick_data(GiNaC::ex Pks_)
      : Pk_string(std::move(Pks_))
      {
      }


  }   // namespace detail
