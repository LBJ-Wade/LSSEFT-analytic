//
// Created by David Seery on 01/09/2017.
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

#include "one_loop_reduced_integral.h"

#include "detail/special_functions.h"
#include "detail/legendre_utils.h"

#include "shared/exceptions.h"
#include "localizations/messages.h"


namespace one_loop_reduced_integral_impl
  {

    integration_element::integration_element(GiNaC::ex ig_, GiNaC::ex ms_, GiNaC::ex wp_, time_function tm_,
                                             GiNaC_symbol_set vs_, GiNaC_symbol_set em_)
      : integrand(std::move(ig_)),
        measure(std::move(ms_)),
        WickProduct(std::move(wp_)),
        tm(std::move(tm_)),
        variables(std::move(vs_)),
        external_momenta(std::move(em_))
      {
      }


    void integration_element::write(std::ostream& str) const
      {
        str << "integral";
        for(const auto& sym : this->variables)
          {
            str << " d" << sym;
          }
        str << '\n';

        str << "  time function = " << this->tm << '\n';
        str << "  measure = " << this->measure << '\n';
        str << "  Wick product = " << this->WickProduct << '\n';
        str << "  integrand = " << this->integrand << '\n';
      }


    void integration_element::simplify(const GiNaC::exmap& map)
      {
        this->integrand = this->integrand.subs(map);
        this->measure = this->measure.subs(map);
        this->WickProduct = this->WickProduct.subs(map);
        this->tm = this->tm.subs(map);
      }


    void integration_element::canonicalize_external_momenta()
      {
        for(const auto& sym : this->external_momenta)
          {
            this->integrand = Legendre_to_cosines(this->integrand, sym);
          }
      }


    GiNaC::ex integration_element::get_UV_limit(unsigned int order) const
      {
        // the total contribution from this element is the product of the integrand, the measure, the
        // Wick product and the time function, but we don't want to Taylor expand the Wick product
        auto prod = this->integrand * this->measure * this->tm;

        // the UV limit occurs when the loop momentum is much greater than any of the external
        // momenta. We can achieve the same result by making a series expansion
        // in the limit that all external momenta are small, ie., by making a Taylor expansion in each
        for(const auto& k : this->external_momenta)
          {
            prod = GiNaC::series_to_poly(prod.expand().series(k, order+1));
          }

        return prod * this->WickProduct;
      }


    key::key(const integration_element& elt_)
      : elt(elt_)
      {
      }


    size_t key::hash() const
      {
        // concatenate symbols in integration variable list
        std::string symbol_string;
        std::for_each(this->elt.variables.begin(), this->elt.variables.end(),
                      [&](const GiNaC::symbol& e) -> std::string
                        { return symbol_string += e.get_name(); });

        // combine both hashes together
        size_t h = 0;
        hash_impl::hash_combine(h, symbol_string);

        // return final value
        return h;
      }


    bool key::is_equal(const key& obj) const
      {
        return std::equal(this->elt.variables.cbegin(), this->elt.variables.cend(),
                          obj.elt.variables.cbegin(), obj.elt.variables.cend(),
                          [](const GiNaC::symbol& asym, const GiNaC::symbol& bsym) -> bool
                            { return asym.get_name() == bsym.get_name(); });
      }


    std::ostream& operator<<(std::ostream& str, const integration_element& obj)
      {
        obj.write(str);
        return str;
      }

  }   // namespace one_loop_reduced_integral_impl


one_loop_reduced_integral::one_loop_reduced_integral(const loop_integral& i_, symbol_factory& sf_)
  : loop_int(i_),
    Rayleigh_momenta(i_.get_Rayleigh_momenta()),
    WickProduct(i_.get_Wick_product()),
    tm(i_.get_time_function()),
    external_momenta(i_.get_external_momenta()),
    sf(sf_),
    x(sf_.make_symbol("x"))
  {
    // throw if we were given a tree-level expression
    if(loop_int.get_loop_order() == 0)
      throw exception(ERROR_ONE_LOOP_REDUCE_WITH_TREE, exception_code::loop_transformation_error);

    // throw if we were given a 2+ loop expression
    if(loop_int.get_loop_order() > 1)
      throw exception(ERROR_ONE_LOOP_REDUCE_WITH_MULTIPLE_LOOPS, exception_code::loop_transformation_error);

    // cache loop momentum
    loop_q = *loop_int.get_loop_momenta().begin();

    // extract momentum kernel from integral
    auto K = loop_int.get_kernel();

    // convert explicit dot products to Cos(a,b) format and expand to get a representation
    // suitable for term-by-term decomposition into a sum of products of Legendre polynomials
    K = dot_products_to_cos(K).expand();

    // apply term-by-term decomposition to K
    if(GiNaC::is_a<GiNaC::mul>(K))
      {
        // just a single product at the top level, so reduce that
        this->reduce(K);
      }
    else if(GiNaC::is_a<GiNaC::add>(K))
      {
        // a sum of terms at top level; work through each one and reduce them term-by-term
        for(size_t i = 0; i < K.nops(); ++i)
          {
            this->reduce(K.op(i));
          }
      }
    else
      throw exception(ERROR_BADLY_FORMED_TOP_LEVEL_MOMENTUM_KERNEL, exception_code::loop_transformation_error);
  }


void one_loop_reduced_integral::reduce(const GiNaC::ex& term)
  {
    // find which Rayleigh momenta this term depends on, if any
    // need to remember that Rayleigh momenta can occur in the momentum kernel but also in the Wick product
    GiNaC_symbol_set Rayleigh_mma;

    for(const auto& rule : this->Rayleigh_momenta)
      {
        const auto& sym = GiNaC::ex_to<GiNaC::symbol>(rule.first);
        if(term.has(sym) || this->WickProduct.has(sym))
          {
            Rayleigh_mma.insert(sym);
          }
      }

    if(Rayleigh_mma.empty())
      {
        // if no Rayleigh momenta then we need only account for the loop momentum
        this->one_loop_reduce_zero_Rayleigh(term);
      }
    else if(Rayleigh_mma.size() == 1)
      {
        // if exactly one Rayleigh momentum then we know how to reduce it
        this->one_loop_reduce_one_Rayleigh(term, *Rayleigh_mma.begin());
      }
    else
      // currently we don't know how to reduce terms with more than one Rayleigh momentum
      throw exception(ERROR_MULTIPLE_RAYLEIGH_MOMENTA_NOT_IMPLEMENTED, exception_code::loop_transformation_error);
  }


void one_loop_reduced_integral::one_loop_reduce_zero_Rayleigh(const GiNaC::ex& term)
  {
    // first, convert all angular terms involving the loop momentum to Legendre representation
    auto temp = Legendre_to_cosines(term, this->loop_q);
    temp = cosines_to_Legendre(temp, this->loop_q);

    // step through expression, identifying terms with zero, one, two or more Legendre polynomials
    // and using the generalized orthogonality relation to perform the angular integrations
    temp = this->integrate_Legendre(temp, this->loop_q,
      [&](auto temp, auto q) -> auto
        {
          return this->apply_Legendre_orthogonality(temp, q);
        });

    // store result if it is nonzero
    if(temp != 0)
      {
        // construct integration element
        auto measure = this->loop_q*this->loop_q / GiNaC::pow(2*GiNaC::Pi, 3);

        using one_loop_reduced_integral_impl::integration_element;
        using one_loop_reduced_integral_impl::key;
        auto elt =
          std::make_unique<integration_element>(temp, measure, this->WickProduct, this->tm,
                                                GiNaC_symbol_set{this->loop_q}, this->external_momenta);

        // insert in database
        this->integrand[key{*elt}].push_back(std::move(elt));
      }
  }


void one_loop_reduced_integral::one_loop_reduce_one_Rayleigh(const GiNaC::ex& term, const GiNaC::symbol& R)
  {
    // we need to perform the integration in an order for which we get only products of two Legendre
    // polynomials involving the integration variable

    // this means we can't do the \hat{x} integral first, because typically there will be three P_n involving \hat{x}.
    // also, we can't guarantee to do the \hat{q} (= loop) integral first because this might involve dot
    // products with \hat{k}, \hat{r} or \hat{x}
    // so we must do the \hat{s} integral first

    // get expansion of the Rayleigh momenta: s = R
    const GiNaC::ex& Rayleigh = this->Rayleigh_momenta.at(R);

    // we represent the \delta function \delta(s-R) in the form exp ix.(s-R), so we need to get the
    // expansion of in terms of the loop momentum and the external momentum
    const GiNaC::ex Rexp = Rayleigh.expand();
    if(Rexp.degree(this->loop_q) == 0)
      throw exception(ERROR_NO_LOOPQ_IN_RAYLEIGH, exception_code::loop_transformation_error);
    if(Rexp.degree(this->loop_q) > 1)
      throw exception(ERROR_DEGREE_OF_LOOPQ_IN_RAYLEIGH_TOO_LARGE, exception_code::loop_transformation_error);
    const auto& loop_coeff = GiNaC::ex_to<GiNaC::numeric>(Rexp.coeff(this->loop_q, 1));
    if(loop_coeff.to_int() != 1 && loop_coeff.to_int() != -1)
      throw exception(ERROR_LOOPQ_HAS_WRONG_COEFF_IN_RAYLEIGH, exception_code::loop_transformation_error);

    using kext_coeff_map = std::map< GiNaC::symbol, GiNaC::numeric >;
    kext_coeff_map kext_coeff;

    for(const auto& k : this->external_momenta)
      {
        if(Rexp.degree(k) == 1)
          {
            const auto& coeff = GiNaC::ex_to<GiNaC::numeric>(Rexp.coeff(k, 1));
            if(coeff.to_int() != 1 && coeff.to_int() != -1)
              throw exception(ERROR_KEXT_HAS_WRONG_COEFF_IN_RAYLEIGH, exception_code::loop_transformation_error);
            kext_coeff[k] = coeff;
          }
        else if(Rexp.degree(k) > 1)
          throw exception(ERROR_DEGREE_OF_KEXT_IN_RAYLEIGH_TOO_LARGE, exception_code::loop_transformation_error);
      }

    if(kext_coeff.empty())
      throw exception(ERROR_NO_KEXT_IN_RAYLEIGH, exception_code::loop_transformation_error);
    if(kext_coeff.size() > 1)
      throw exception(ERROR_TOO_MANY_KEXT_IN_RAYLEIGH, exception_code::loop_transformation_error);

    // first task is to perform the Rayleigh angular integration
    // there should not be any angular dependence in the momentum kernel, making the integral trivial
    auto temp = Legendre_to_cosines(term, R);
    temp = cosines_to_Legendre(temp, R);

    auto R_pairs = get_LegP_pairs(temp, R);
    if(!R_pairs.empty())
      throw exception(ERROR_KERNEL_DEPENDS_ON_ANGULAR_RAYLEIGH_MOMENTUM, exception_code::loop_transformation_error);

    // since the integral is trivial we just get a factor of 4pi; factors of (2l+1) are all unity
    temp *= GiNaC::numeric{4} * GiNaC::Pi;

    // at this stage we can do the \hat{x} integral
    // it will couple together the angular terms in the Rayleigh plane wave expansion for k and L
    // (with suitable factors of 4pi/(2l+1)
    // this doesn't make any difference to the kernel, as long as we remember that it has been done

    // it generates a sum of the form
    // sum_{n=1}^\infty (-1)^n 4pi (2n+1) LegP(n, k.L)
    // times a sign factor generated by the sign of k and L

    temp = Legendre_to_cosines(temp, this->loop_q);
    temp = cosines_to_Legendre(temp, this->loop_q);

    // step through this expression, pairing up LegP(n, k.L) and LegP(n, r.L) terms,
    // either appearing explicitly or from the sum generated by the \hat{x} integral,
    // and also replacing the dx integral with a Fabrikant function
    temp = this->integrate_Legendre(temp, this->loop_q,
      [&](auto term, auto q) -> auto
        {
          return this->apply_Legendre_orthogonality(term, q, loop_coeff, kext_coeff.begin()->first, kext_coeff.begin()->second, R);
        });

    // the input kernel should be dimensionless
    // the output kernel has an extra integral d^3 s and should therefore have dimension [k^-3], ie. it loses
    // three powers of momentum compared to the input kernel

    // store result if it is nonzero
    if(temp != 0)
      {
        // construct integration element, assuming that the Fabrikant integrals imposed the
        // triangle condition on R, loop_q and kext
        auto measure = this->loop_q*this->loop_q / GiNaC::pow(2*GiNaC::Pi, 3);

        auto& kext = kext_coeff.begin()->first;
        GiNaC::ex R_replace = GiNaC::sqrt(kext*kext + this->loop_q*this->loop_q - 2*this->loop_q*kext*this->x);

        // measure is R^2 dR, and dR -> k q / R
        measure *= R_replace * kext * this->loop_q / GiNaC::pow(2*GiNaC::Pi, 3);

        GiNaC::exmap R_map;
        R_map[R] = R_replace;
        temp = temp.subs(R_map);

        using one_loop_reduced_integral_impl::integration_element;
        using one_loop_reduced_integral_impl::key;
        auto elt =
          std::make_unique<integration_element>(temp, measure, this->WickProduct, this->tm,
                                                GiNaC_symbol_set{this->loop_q, this->x}, this->external_momenta);

        // insert in database
        this->integrand[key{*elt}].push_back(std::move(elt));
      }
  }


GiNaC::ex one_loop_reduced_integral::apply_Legendre_orthogonality(const GiNaC::ex& expr, const GiNaC::symbol& q)
  {
    GiNaC::ex temp{1};

    if(!GiNaC::is_a<GiNaC::mul>(expr))
      throw exception(ERROR_BADLY_FORMED_LEGENDRE_SUM_TERM, exception_code::loop_transformation_error);

    using Legendre_list = std::vector< std::pair< GiNaC::symbol, unsigned int > >;
    Legendre_list partner_q;

    // run through each factor in the expression
    // if it is a Legendre polynomial involving q then record the momentum it occurs with, and its order
    for(size_t i = 0; i < expr.nops(); ++i)
      {
        const auto term = expr.op(i);
        if(!GiNaC::is_a<GiNaC::function>(term)) { temp *= term; continue; }

        const auto& fn = GiNaC::ex_to<GiNaC::function>(term);
        if(fn.get_name() != "LegP") { temp *= term; continue; }

        auto n = static_cast<unsigned int>(GiNaC::ex_to<GiNaC::numeric>(fn.op(0)).to_int());
        auto p1 = GiNaC::ex_to<GiNaC::symbol>(fn.op(1));
        auto p2 = GiNaC::ex_to<GiNaC::symbol>(fn.op(2));

        if(p1 == q) { partner_q.emplace_back(p2, n); continue; }
        if(p2 == q) { partner_q.emplace_back(p1, n); continue; }

        // not a Legendre polynomial involving q, so move on
        temp *= term;
      }

    // if no Legendre polynomials, equivalent to LegP(0, x)
    if(partner_q.empty())
      {
        return GiNaC::numeric{4} * GiNaC::Pi * temp;
      }

    // if one Legendre polynomial, gives nonzero if n=0
    if(partner_q.size() == 1)
      {
        if(partner_q.front().second != 0) return GiNaC::ex{0};
        return GiNaC::numeric{4} * GiNaC::Pi * temp;
      }

    // if more than two Legendre polynomials, don't know what to do
    if(partner_q.size() > 2)
        throw exception(ERROR_CANT_INTEGRATE_MORE_THAN_TWO_LEGP, exception_code::loop_transformation_error);

    // remaining case is two Legendre polynomials, which can be integrated
    // using the orthogonality relation
    if(partner_q.front().second != partner_q.back().second) return GiNaC::ex{0};

    unsigned int n = partner_q.front().second;
    auto p1 = partner_q.front().first;
    auto p2 = partner_q.back().first;

    if(std::less<GiNaC::symbol>{}(p2, p1)) std::swap(p1, p2);

    auto cf = GiNaC::numeric{4} * GiNaC::Pi / (2*GiNaC::numeric{n} + 1);
    return cf * (n > 0 ? Angular::LegP(n, p1, p2).expand() : GiNaC::numeric{1}.expand()) * temp;
  }


static std::map< unsigned int, GiNaC::numeric > A_cache;


GiNaC::numeric A(unsigned int r)
  {
    if(A_cache.find(r) != A_cache.end()) return A_cache[r];

    unsigned int numerator = 1;
    unsigned int count = 2*r - 1;

    while(count > 1)
      {
        numerator *= count;
        count -= 2;
      }

    unsigned int denominator = 1;
    count = r;

    while(count > 1)
      {
        denominator *= count;
        count -= 1;
      }

    A_cache[r] = GiNaC::numeric{numerator} / GiNaC::numeric{denominator};
    return A_cache[r];
  }


// compute the integral over L of
// [ sum_n (-1)^n 4pi (2n+1) LegP(n, k.L) ] LegP[p, k.L] LegP[q, r.L]
// using the Neumann-Adams formula to handle the product of Legendre polynomials of k.L
GiNaC::ex NeumannAdamsSum(const GiNaC::symbol& L, const GiNaC::numeric& Lcoeff, const GiNaC::symbol& k,
                          const GiNaC::numeric& kcoeff, unsigned int p, const GiNaC::symbol& r, unsigned int q,
                          const GiNaC::symbol& R)
  {
    GiNaC::ex expr{0};

    // Neumann-Adams product of LegP(n, k.L) and LegP(p, k.L) will generate terms between n+p and |n-p|.
    // These terms have to match LegP(q, r.L) which requires |n-p| <= q <= n+p, or
    // n-p <= q --> n <= p+q
    // p-n <= q --> n >= p-q
    // q <= n+p --> n >= q-p
    // so |p-q| <= n <= p+q

    auto n_lo = static_cast<unsigned int>(std::abs(static_cast<int>(p) - static_cast<int>(q)));
    auto n_hi = p + q;

    for(unsigned int n = n_lo; n <= n_hi; ++n)
      {
        // determine whether there is a value of r that will make n + p - 2t = q, with 0 <= r <= min(n,p)
        auto twot = static_cast<int>(n) + static_cast<int>(p) - static_cast<int>(q);

        // can't solve this equation if 2r is not even
        if(twot % 2 != 0) continue;
        auto t = twot / 2;

        if(t < 0) continue;
        if(t > std::min(n, p)) continue;

        auto a = A(n - t) * A(t) * A(p - t) / A(n + p - t);
        auto b = GiNaC::numeric{2*static_cast<int>(n) + 2*static_cast<int>(p) - 4*t + 1};
        auto c = GiNaC::numeric{2*static_cast<int>(n) + 2*static_cast<int>(p) - 2*t + 1};
        auto coeff = a * b / c;

        // in cfn two factors of (2l+1) in numerator from Rayleigh expansion compete with (2l+1) in denominator
        // to give (2l+1) in numerator
        // cfq has contributions only from orthogonality relation
        auto cfq = GiNaC::numeric{4} * GiNaC::Pi / (2*GiNaC::numeric{q} + 1);
        auto cfn = GiNaC::numeric{4} * GiNaC::Pi * (2*GiNaC::numeric{n} + 1);

        auto p1 = k;
        auto p2 = r;
        if(std::less<GiNaC::symbol>{}(p2, p1)) std::swap(p1, p2);

        // to keep track of the k and L coefficients, remember that these are attached to the sum from the
        // Rayleigh expansion, so these can produce a sign flip *only* with the term Leg(n, k.L)
        expr += cfq*cfn * coeff * GiNaC::pow(Lcoeff, n) * GiNaC::pow(kcoeff, n) * GiNaC::pow(-1, n)
                * (q > 0 ? Angular::LegP(q, p1, p2).expand() : GiNaC::numeric{1}.expand()) * Fabrikant::FabJ(0, n, n, R, L, k);
      }

    return expr;
  }

GiNaC::ex
one_loop_reduced_integral::apply_Legendre_orthogonality(const GiNaC::ex& expr, const GiNaC::symbol& L, const GiNaC::numeric& Lcoeff,
                                                        const GiNaC::symbol& k, const GiNaC::numeric& kcoeff, const GiNaC::symbol& R)
  {
    GiNaC::ex temp{1};

    using Legendre_list = std::vector< std::pair< GiNaC::symbol, unsigned int > >;
    Legendre_list partner_q;

    for(size_t i = 0; i < expr.nops(); ++i)
      {
        const auto term = expr.op(i);
        if(!GiNaC::is_a<GiNaC::function>(term)) { temp *= term; continue; }

        const auto& fn = GiNaC::ex_to<GiNaC::function>(term);
        if(fn.get_name() != "LegP") { temp *= term; continue; }

        auto n = static_cast<unsigned int>(GiNaC::ex_to<GiNaC::numeric>(fn.op(0)).to_int());
        auto p1 = GiNaC::ex_to<GiNaC::symbol>(fn.op(1));
        auto p2 = GiNaC::ex_to<GiNaC::symbol>(fn.op(2));

        if(p1 == L) { partner_q.emplace_back(p2, n); continue; }
        if(p2 == L) { partner_q.emplace_back(p1, n); continue; }

        // not a Legendre polynomial involving L, so move on
        temp *= term;
      }

    // if no Legendre polynomials, only the LegP(0, k.L) term in the sum contributes
    if(partner_q.empty()) return NeumannAdamsSum(L, Lcoeff, k, kcoeff, 0, k, 0, R) * temp;

    // if one Legendre polynomial, depends whether the partner field is k or something else
    if(partner_q.size() == 1)
      {
        if(partner_q.front().first == k)
          {
            // this is another Legendre polynomial of the form LegP(m, k.L), so we need to expand the product
            // LegP(n, k.L) * LegP(n, k.L) and integrate against 1 = Leg(0, r.L)
            return NeumannAdamsSum(L, Lcoeff, k, kcoeff, partner_q.front().second, k, 0, R) * temp;
          }
        else
          {
            // this is a Legendre polynomial of the form LegP(m, r.L). This is a special case of the
            // Neumann-Adams sum too
            return NeumannAdamsSum(L, Lcoeff, k, kcoeff, 0, partner_q.front().first, partner_q.front().second, R) * temp;
          }
      }

    // if more than two Legendre polynomials, don't know what to do
    if(partner_q.size() > 2)
      throw exception(ERROR_CANT_INTEGRATE_MORE_THAN_TWO_LEGP, exception_code::loop_transformation_error);

    if(static_cast<bool>(partner_q.front().first == k) && static_cast<bool>(partner_q.back().first == k))
      throw exception(ERROR_TOO_MANY_KEXT_IN_LEGENDRE_SUM, exception_code::loop_transformation_error);

    if(partner_q.front().first == k)
      return NeumannAdamsSum(L, Lcoeff, k, kcoeff, partner_q.front().second, partner_q.back().first, partner_q.back().second, R) * temp;

    if(partner_q.back().first == k)
      return NeumannAdamsSum(L, Lcoeff, k, kcoeff, partner_q.back().second, partner_q.front().first, partner_q.front().second, R) * temp;

    throw exception(ERROR_FAILED_TO_REDUCE_TRIPLE_PRODUCT_LEGP, exception_code::loop_transformation_error);
  }


void one_loop_reduced_integral::write(std::ostream& out) const
  {
    using one_loop_reduced_integral_impl::integration_element;

    for(const auto& record : this->integrand)
      {
        const auto& vars = record.first;
        const auto& data = record.second;

        for(const auto& iptr : data)
          {
            const auto& i = *iptr;
            out << i;
          }
      }
  }


void one_loop_reduced_integral::simplify(const GiNaC::exmap& map)
  {
    using one_loop_reduced_integral_impl::integration_element;

    for(const auto& record : this->integrand)
      {
        const auto& data = record.second;

        for(const auto& iptr : data)
          {
            auto& i = *iptr;
            i.simplify(map);
          }
      }
  }


void one_loop_reduced_integral::canonicalize_external_momenta()
  {
    using one_loop_reduced_integral_impl::integration_element;

    for(const auto& record : this->integrand)
      {
        const auto& data = record.second;

        for(const auto& iptr : data)
          {
            auto& i = *iptr;
            i.canonicalize_external_momenta();
          }
      }
  }


GiNaC::ex one_loop_reduced_integral::get_UV_limit(unsigned int order) const
  {
    using one_loop_reduced_integral_impl::integration_element;

    GiNaC::ex expr{0};

    for(const auto& record : this->integrand)
      {
        const auto& data = record.second;

        for(const auto& iptr : data)
          {
            auto& i = *iptr;
            auto term = i.get_UV_limit(order);

            const auto& ivars = i.get_integration_variables();

            // if the integration variables for i contain x, then
            // perform a symbolic integration for it
            auto t = ivars.find(this->x);
            if(t != ivars.end())
              {
                term = GiNaC::integral(this->x, -1, 1, term).eval_integ();
              }

            expr += term;
          }
      }

    return expr;
  }


std::ostream& operator<<(std::ostream& str, const one_loop_reduced_integral& obj)
  {
    obj.write(str);
    return str;
  }
