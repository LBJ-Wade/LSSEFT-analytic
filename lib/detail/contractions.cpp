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

    //! Taylor expansion of Pk needs special attention, because we don't want to end up doing
    //! an expansion of something Pk(k).
    //! On the other hand, Pk(sqrt(k^2 + K^2 - 2kLx)) can certainly be expanded
    GiNaC::ex Pk_series(const GiNaC::ex& a, const GiNaC::ex& b, const GiNaC::ex& x,
                        const GiNaC::relational& rel, int order, unsigned options);

    // register Pk as a dummy function
    REGISTER_FUNCTION(Pk, series_func(Pk_series))

    GiNaC::ex Pk_series(const GiNaC::ex& a, const GiNaC::ex& b, const GiNaC::ex& x,
                        const GiNaC::relational& rel, int order, unsigned options)
      {
        // get symbol that we're expanding with respect to
        const auto& sym = GiNaC::ex_to<GiNaC::symbol>(rel.lhs());

        // if series expansion of argument has no constant term, then we don't want to do a Taylor expansion
        if(x.coeff(sym, 0) != 0) throw GiNaC::do_taylor();

        return GiNaC::pseries(rel, GiNaC::epvector{ {Pk(a, b, x), 0} });
      }

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


  }   // namespace detail
