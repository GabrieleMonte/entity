#ifndef PROBLEM_GENERATOR_H
#define PROBLEM_GENERATOR_H

#include "enums.h"
#include "global.h"

#include "arch/kokkos_aliases.h"
#include "traits/pgen.h"
#include "utils/error.h"
#include "utils/numeric.h"

#include "archetypes/particle_injector.h"
#include "archetypes/utils.h"
#include "framework/domain/domain.h"
#include "framework/domain/metadomain.h"
#include "framework/parameters/parameters.h"

#include <string>
#include <utility>
#include <vector>

namespace user {
  using namespace ntt;
  using prmvec_t = std::vector<real_t>;

  /*
   * Dark-photon dark matter resonant conversion, following
   * Hook, Huang & Shalaby, arXiv:2510.13956.
   *
   * The dark photon enters as an external, spatially uniform (k = 0),
   * time-oscillating electric field acting on the Lorentz force only -- it is
   * never a source term in Maxwell's equations. All feedback on the drive is
   * through the deposited current.
   */


  /*
   * ------------------------------------------------------------------
   *  Quiet start (deterministic loading)
   * ------------------------------------------------------------------
   * Entity's stock injector places particles at RANDOM sub-cell positions
   * (src/kernels/injectors.hpp:202). Co-locating the pair already gives
   * rho(x,0) = 0, but the TOTAL density n_e + n_i still fluctuates at
   * ~sqrt(2/ppc), and that is what spreads omega_p(x) and detunes the
   * resonance. The paper loads particles on an even lattice instead.
   */

  /*
    Inverse standard-normal CDF: Acklam's rational approximation (~1.15e-9)
    refined by one Halley step against erfc, which takes it to ~1e-15.
    A 1D Maxwellian of temperature T has u_x ~ Normal(0, sqrt(T)), so
    u = sqrt(T) * probit(q) maps a quantile to a velocity.
  */
  Inline auto probit(real_t q) -> real_t {
    const real_t a0 = -3.969683028665376e+01, a1 = 2.209460984245205e+02,
                 a2 = -2.759285104469687e+02, a3 = 1.383577518672690e+02,
                 a4 = -3.066479806614716e+01, a5 = 2.506628277459239e+00;
    const real_t b0 = -5.447609879822406e+01, b1 = 1.615858368580409e+02,
                 b2 = -1.556989798598866e+02, b3 = 6.680131188771972e+01,
                 b4 = -1.328068155288572e+01;
    const real_t c0 = -7.784894002430293e-03, c1 = -3.223964580411365e-01,
                 c2 = -2.400758277161838e+00, c3 = -2.549732539343734e+00,
                 c4 = 4.374664141464968e+00,  c5 = 2.938163982698783e+00;
    const real_t d0 = 7.784695709041462e-03,  d1 = 3.224671290700398e-01,
                 d2 = 2.445134137142996e+00,  d3 = 3.754408661907416e+00;
    const real_t plo = 0.02425, sqrt2 = 1.4142135623730951,
                 sqrt2pi = 2.5066282746310002;
    real_t x;
    if (q < plo) {
      const real_t t = math::sqrt(-TWO * math::log(q));
      x = (((((c0 * t + c1) * t + c2) * t + c3) * t + c4) * t + c5) /
          ((((d0 * t + d1) * t + d2) * t + d3) * t + ONE);
    } else if (q <= ONE - plo) {
      const real_t t = q - HALF, r = t * t;
      x = (((((a0 * r + a1) * r + a2) * r + a3) * r + a4) * r + a5) * t /
          (((((b0 * r + b1) * r + b2) * r + b3) * r + b4) * r + ONE);
    } else {
      const real_t t = math::sqrt(-TWO * math::log(ONE - q));
      x = -(((((c0 * t + c1) * t + c2) * t + c3) * t + c4) * t + c5) /
           ((((d0 * t + d1) * t + d2) * t + d3) * t + ONE);
    }
    const real_t e = HALF * math::erfc(-x / sqrt2) - q;
    const real_t u = e * sqrt2pi * math::exp(x * x * HALF);
    return x - u / (ONE + x * u * HALF);
  }

  /*
    The stratified set {probit((j+1/2)/ppc)} has mean 0 by symmetry but variance
    slightly BELOW 1: the midpoint rule under-samples the Gaussian tails. The
    deficit is 2.0% at ppc=64 and 0.64% at ppc=200 -- i.e. the plasma would come
    out systematically colder than the temperature asked for, which is exactly
    the sort of silent bias a quiet start exists to avoid. Measure the realized
    variance and divide it out, so the loaded temperature is exact.
  */
  /*
    splitmix64 -> uniform in (0,1). Deterministic and keyed on the GLOBAL cell
    index, so the loading is identical however MPI splits the domain -- the same
    property the lattice loader had, kept here.
  */
  Inline auto hash_u01(std::uint64_t z) -> real_t {
    z += 0x9E3779B97F4A7C15ull;
    z  = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z  = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z  =  z ^ (z >> 31);
    // (z>>11 + 1/2)/2^53 stays strictly inside (0,1) so probit is finite
    return (static_cast<real_t>(z >> 11) + HALF) /
           static_cast<real_t>(9007199254740992.0);   // 2^53
  }

  struct QuietNorm_kernel {
    const npart_t ppc;

    explicit QuietNorm_kernel(npart_t ppc) : ppc { ppc } {}

    Inline void operator()(npart_t j, real_t& acc) const {
      const real_t u = probit((static_cast<real_t>(j) + HALF) /
                              static_cast<real_t>(ppc));
      acc += u * u;
    }
  };

  /*
    Deterministic loading kernel. One thread per electron; the ion is placed at
    the SAME position, so the charge density is 0 bit-for-bit at t=0.

    THE PROPERTY THAT MATTERS is not that t=0 looks clean -- it is that the
    order SURVIVES streaming. Index by velocity CLASS, not by sub-cell slot:

        c   = p / ncell            velocity class, exactly ncell particles each
        i   = p % ncell            cell
        phi = ((c*P0) mod ppc + 1/2)/ppc          depends ONLY on c
        x   = xi_min + i + phi

    so every velocity class is a perfect lattice of spacing dx spanning the box.
    Under free streaming a class only translates, and a dx-spaced lattice
    deposits exactly uniform density at ANY translation (the shape function is a
    partition of unity at that spacing) -- so the total density stays uniform to
    round-off rather than decaying to the random floor.

    Tying the velocity to the sub-cell slot instead (with a per-cell rotation)
    fails exactly here: a given velocity then sits at a different offset in each
    cell, the classes shear apart as soon as particles move, and the noise
    reappears within one step. Measured that way: dn/n = 6.3e-3 after a single
    step, and a floor only 2x below the random start.

    The coprime stride P0 scrambles which offset gets which velocity, so no cell
    carries a monotonic v(x) ramp. Ions reuse the positions (co-location is what
    makes the charge density vanish) but take their velocities through a second
    stride, so the two species are not phase-space correlated.

    Nothing here depends on the cell index, so the loading is identical however
    MPI decomposes the domain: an 8-node run reproduces a 1-node run.
  */
  template <SimEngine::type S, class M>
  struct QuietInit_kernel {
    static constexpr auto D = M::Dim;

    array_t<int*>      i1_1, i1_2;
    array_t<prtldx_t*> dx1_1, dx1_2;
    array_t<real_t*>   ux1_1, ux2_1, ux3_1, w_1;
    array_t<real_t*>   ux1_2, ux2_2, ux3_2, w_2;
    array_t<short*>    tag_1, tag_2;

    const array_t<real_t*> xi_min;
    const npart_t          ppc, ncell;
    const std::size_t      P0, P2;
    const real_t           sig1, sig2;
    const bool             one_v;
    const bool             shear;

    QuietInit_kernel(Particles<D, M::CoordType>& sp1,
                     Particles<D, M::CoordType>& sp2,
                     const array_t<real_t*>& xi_min, npart_t ppc, npart_t ncell,
                     std::size_t P0, std::size_t P2,
                     real_t sig1, real_t sig2, bool one_v, bool shear)
      : i1_1 { sp1.i1 }, i1_2 { sp2.i1 }
      , dx1_1 { sp1.dx1 }, dx1_2 { sp2.dx1 }
      , ux1_1 { sp1.ux1 }, ux2_1 { sp1.ux2 }, ux3_1 { sp1.ux3 }, w_1 { sp1.weight }
      , ux1_2 { sp2.ux1 }, ux2_2 { sp2.ux2 }, ux3_2 { sp2.ux3 }, w_2 { sp2.weight }
      , tag_1 { sp1.tag }, tag_2 { sp2.tag }
      , xi_min { xi_min }, ppc { ppc }, ncell { ncell }
      , P0 { P0 }, P2 { P2 }, sig1 { sig1 }, sig2 { sig2 }, one_v { one_v }
      , shear { shear } {}

    Inline void operator()(npart_t p) const {
      real_t  phi, q1, q2;
      npart_t i;

      if (shear) {
        /*
          POSITIONS: an exact lattice. Cell i holds ppc particles at sub-cell
          offsets (s+1/2)/ppc, so density is uniform to round-off and
          eps_noise(t=0) = 0 -- the paper's stated initial condition.

          VELOCITIES: pseudo-random per particle. The paper constrains this only
          implicitly ("all particles equally spaced, with ions and electrons
          colocated" describes POSITIONS), and it is what lets the lattice
          decohere into physical shot noise, as their Fig. A00 shows happening
          within a few 1/omega_p.

          A rigid per-cell rotation of a stratified set does NOT work here: each
          velocity class would still be a uniform comb, of spacing
          (1 - R/ppc)*dx rather than dx, and a uniform comb translates rigidly
          and never clumps. The assignment has to be irregular cell to cell.
        */
        i = p / ppc;
        const npart_t s = p % ppc;
        phi = (static_cast<real_t>(s) + HALF) / static_cast<real_t>(ppc);
        const std::uint64_t gcell =
          static_cast<std::uint64_t>(xi_min(0) + static_cast<real_t>(i));
        const std::uint64_t key = gcell * 1000003ull + static_cast<std::uint64_t>(s);
        q1 = hash_u01(key * 4ull);
        q2 = hash_u01(key * 4ull + 1ull);
      } else {
        /*
          LEGACY "quiet_lattice": velocity keyed on the CLASS, so every class is
          a perfect dx-spaced lattice that survives streaming forever. Kept only
          to reproduce the earlier runs. It does not decohere, seeds daughter
          modes from round-off instead of 1/sqrt(N_pc), and therefore does NOT
          reproduce the paper.
        */
        const npart_t     c  = p / ncell;
        i = p % ncell;
        const std::size_t np = static_cast<std::size_t>(ppc);
        const std::size_t cc = static_cast<std::size_t>(c);
        phi = (static_cast<real_t>((cc * P0) % np) + HALF) /
              static_cast<real_t>(ppc);
        const real_t inv = ONE / static_cast<real_t>(ppc);
        q1 = (static_cast<real_t>(c) + HALF) * inv;
        q2 = (static_cast<real_t>((cc * P2) % np) + HALF) * inv;
      }

      const real_t   x_Cd = xi_min(0) + static_cast<real_t>(i) + phi;
      const int      i1   = static_cast<int>(x_Cd);
      const prtldx_t dx1  = static_cast<prtldx_t>(x_Cd - static_cast<real_t>(i1));

      i1_1(p) = i1;  dx1_1(p) = dx1;
      ux1_1(p) = sig1 * probit(q1);
      ux2_1(p) = ZERO;  ux3_1(p) = ZERO;
      w_1(p) = ONE;  tag_1(p) = ParticleTag::alive;

      i1_2(p) = i1;  dx1_2(p) = dx1;   // co-located: charge density 0 exactly
      ux1_2(p) = sig2 * probit(q2);
      ux2_2(p) = ZERO;  ux3_2(p) = ZERO;
      w_2(p) = ONE;  tag_2(p) = ParticleTag::alive;

      if (not one_v) {
        // 3V transverse components; they feel no force in a 1D electrostatic
        // run, so these affect only the T00 bookkeeping.
        const std::size_t np  = static_cast<std::size_t>(ppc);
        const std::size_t cc  = static_cast<std::size_t>(p / ncell);
        const real_t      inv = ONE / static_cast<real_t>(ppc);
        const std::uint64_t gcell =
          static_cast<std::uint64_t>(xi_min(0) + static_cast<real_t>(p / ppc));
        const std::uint64_t key =
          gcell * 1000003ull + static_cast<std::uint64_t>(p % ppc);
        const real_t q1y = shear
          ? hash_u01(key * 4ull + 2ull)
          : (static_cast<real_t>((cc * P0 + 5u) % np) + HALF) * inv;
        const real_t q1z = shear
          ? hash_u01(key * 4ull + 3ull)
          : (static_cast<real_t>((cc * P2 + 7u) % np) + HALF) * inv;
        ux2_1(p) = sig1 * probit(q1y);
        ux3_1(p) = sig1 * probit(q1z);
        ux2_2(p) = sig2 * probit(q1z);
        ux3_2(p) = sig2 * probit(q1y);
      }
    }
  };

  template <Dimension D>
  struct ExtFields {

    /*
      @param e1: the (uniform) external electric field in units of B0
    */
    ExtFields(real_t e1) : e1 { e1 } {}

    // e_ext: external electric field, uniform in space
    Inline auto ex1(const coord_t<D>&) const -> real_t {
      return e1;
    }

  private:
    const real_t e1;
  };

  template <SimEngine::type S, class M>
  struct PGen {
    static constexpr auto D { M::Dim };
    // compatibility traits for the problem generator
    static constexpr auto engines = ::traits::pgen::compatible_with<SimEngine::SRPIC> {};
    static constexpr auto metrics =
      ::traits::pgen::compatible_with<Metric::Minkowski> {};
    static constexpr auto dimensions =
      ::traits::pgen::compatible_with<Dim::_1D, Dim::_2D, Dim::_3D> {};

    const SimulationParams& params;

    // drive parameters
    const std::string drive, loading;
    const bool        one_v;
    const std::string sweep_phase;
    const real_t      A0, omega, omega0, domega_dt;
    // plasma parameters
    const prmvec_t    densities, temperatures;

    PGen(const SimulationParams& p, const Metadomain<S, M>& global_domain)
      : params { p }
      , drive { p.template get<std::string>("setup.drive", "resonance") }
      , loading { p.template get<std::string>("setup.loading", "random") }
      // 1V matches SHARP exactly: the transverse DOF feel no force in a 1D
      // electrostatic run, so zeroing them changes only the T00 bookkeeping
      // (1/2 kT per species rather than 3/2 kT) and removes the 3V-vs-1V
      // correction from every downstream comparison.
      , one_v { p.template get<bool>("setup.one_v", true) }
      // See drive_phase() below. "integrated" is correct and is the default;
      // "naive" exists to test the literal reading of the paper's notation.
      , sweep_phase { p.template get<std::string>("setup.sweep_phase",
                                                 "integrated") }
      , A0 { p.template get<real_t>("setup.A0", ZERO) }
      , omega { p.template get<real_t>("setup.omega", ONE) }
      , omega0 { p.template get<real_t>("setup.omega0", ZERO) }
      , domega_dt { p.template get<real_t>("setup.domega_dt", ZERO) }
      , densities { p.template get<prmvec_t>("setup.densities", prmvec_t { TWO }) }
      , temperatures { p.template get<prmvec_t>("setup.temperatures",
                                                prmvec_t { ZERO, ZERO }) } {
      const auto nspec = params.template get<std::size_t>("particles.nspec");
      raise::ErrorIf(nspec != 2,
                     "This setup requires exactly two species (electrons + ions)",
                     HERE);
      raise::ErrorIf(
        global_domain.species_params()[0].charge() !=
          -global_domain.species_params()[1].charge(),
        "Charges of the two species must be opposite for this setup",
        HERE);
      raise::ErrorIf(drive != "resonance" and drive != "landau_zener",
                     "setup.drive must be either `resonance` or `landau_zener`",
                     HERE);
      raise::ErrorIf(loading != "random" and loading != "quiet" and
                       loading != "quiet_lattice",
                     "setup.loading must be `random`, `quiet`, or `quiet_lattice`",
                     HERE);
      raise::ErrorIf(sweep_phase != "integrated" and sweep_phase != "naive",
                     "setup.sweep_phase must be either `integrated` or `naive`",
                     HERE);
      raise::ErrorIf(densities.size() != 1,
                     "setup.densities must have a single entry (per species pair)",
                     HERE);
      raise::ErrorIf(temperatures.size() != 2,
                     "setup.temperatures must have one entry per species",
                     HERE);
    }

    /*
      Two conventions, selected by setup.sweep_phase.

      "integrated" (default, and the correct one): phi = int_0^t omega(t') dt'
          = omega0*t + (1/2)*domega_dt*t^2.  The instantaneous frequency
          d(phi)/dt is then omega0 + domega_dt*t, which is what "a drive whose
          frequency evolves as omega(t)" means. With omega0 = 0.8 and
          domega_dt = 2e-5 the crossing omega = omega_p is at omega_p t = 1e4.

      "naive": phi = omega(t)*t = omega0*t + domega_dt*t^2, the literal reading
          of the paper's "A_0 cos(omega t)" with "omega(t) = 0.8 w_p +
          0.1 w_p t/5000".  This is a well-defined chirp but its instantaneous
          frequency is omega0 + 2*domega_dt*t -- TWICE the nominal sweep rate --
          so the crossing moves to omega_p t = 5000.

      The paper's internal evidence is contradictory. Its stated omega_p/H ~ 7e4
      matches "integrated": with d ln(omega_p^2)/dt = 3H and H = 2*domega_dt/(3*omega)
      at crossing, domega_dt = 2e-5 gives omega_p/H = 7.5e4, whereas the naive
      form's effective 4e-5 gives 3.75e4. But its prose places the resonance at
      omega_p t = 5000, which is where "naive" crosses. (That prose is weak
      evidence: LZ transfer peaks well before the crossing anyway -- our own
      integrated run peaks at omega_p t = 3689 with the crossing at 1e4.)

      "naive" is provided to measure how sensitive the nonlinear heating is to
      the sweep rate, not because it is right.
    */
    auto drive_phase(simtime_t time) const -> real_t {
      if (drive == "landau_zener") {
        if (sweep_phase == "naive") {
          return (omega0 + domega_dt * time) * time;
        }
        return omega0 * time + HALF * domega_dt * time * time;
      } else {
        return omega * time;
      }
    }

    auto ext_efield(simtime_t time) const -> real_t {
      return A0 * math::cos(drive_phase(time));
    }

    /*
      @returns a pair of (apply_external_fields, external_fields)

      @note applied to both species: the pusher scales by (q/m) per species, so
            the ions automatically receive 1/1836 of the electron kick
      @note the drive is uniform in space, so cos() is evaluated once per
            species per step here rather than once per particle in the kernel
    */
    auto ExternalFields(simtime_t time, spidx_t, const Domain<S, M>&) const
      -> std::pair<bool, ExtFields<M::Dim>> {
      return {
        A0 != ZERO,
        ExtFields<M::Dim> { ext_efield(time) }
      };
    }

    /*
      Records the instantaneous drive in the stats output, so that the amplitude
      and phase of the drive can be checked against the response directly.
    */
    auto CustomStat(const std::string&  name,
                    timestep_t,
                    simtime_t           time,
                    const Domain<S, M>&) const -> real_t {
      if (name == "e_ext") {
        return ext_efield(time);
      } else {
        raise::Error("Unrecognized custom stat: " + name, HERE);
        return ZERO;
      }
    }

    /*
      Uniform Maxwellian electrons and ions, injected as a co-located pair, so
      that rho(x, t = 0) = 0 exactly.
    */
    void InitPrtls(Domain<S, M>& domain) {
      if (loading == "random") {
        arch::InjectUniformMaxwellians<S, M>(params,
                                             domain,
                                             densities[0],
                                             { temperatures[0], temperatures[1] },
                                             { 1, 2 });
        return;
      }

      // ---------------- quiet start ----------------
      boundaries_t<real_t> box;
      for (auto d { 0u }; d < M::Dim; ++d) {
        box.push_back(Range::All);
      }
      const auto result = arch::ComputeNumInject<S, M>(params, domain,
                                                       densities[0], box);
      if (not std::get<0>(result)) {
        return;
      }
      const auto npart  = std::get<1>(result);
      const auto xi_min = std::get<2>(result);
      const auto xi_max = std::get<3>(result);

      auto xi_min_h = Kokkos::create_mirror_view(xi_min);
      auto xi_max_h = Kokkos::create_mirror_view(xi_max);
      Kokkos::deep_copy(xi_min_h, xi_min);
      Kokkos::deep_copy(xi_max_h, xi_max);

      const auto ncell = static_cast<npart_t>(xi_max_h(0) - xi_min_h(0) + HALF);
      raise::ErrorIf(ncell == 0u, "Quiet start: empty local domain", HERE);
      raise::ErrorIf(npart % ncell != 0u,
                     "Quiet start: particles/cell must be an integer; check "
                     "that ppc0 * densities[0] * 0.5 is a whole number",
                     HERE);
      const auto ppc = npart / ncell;
      raise::ErrorIf(ppc < 2u, "Quiet start: need at least 2 particles/cell", HERE);

      // stride coprime to ppc, near ppc/phi: scatters consecutive slots across
      // the quantile range so no cell carries a systematic v(x) ramp
      // BOTH strides must be coprime to ppc itself -- that is what makes
      // j -> (j*P + shift) mod ppc a bijection, hence what guarantees each cell
      // carries the complete stratified quantile set exactly once.
      const auto coprime_near = [](std::size_t n, long double frac) -> std::size_t {
        auto gcd = [](std::size_t a, std::size_t b) {
          while (b != 0u) { const auto t = b; b = a % b; a = t; }
          return a;
        };
        auto k = static_cast<std::size_t>(static_cast<long double>(n) * frac);
        if (k < 2u) { k = 2u; }
        while (k > 1u and gcd(k, n) != 1u) { --k; }
        return (k > 1u) ? k : 1u;
      };
      const auto np = static_cast<std::size_t>(ppc);
      const auto P0 = coprime_near(np, 0.6180339887498949L);   // 1/phi
      auto       P2 = coprime_near(np, 0.3819660112501051L);   // 1/phi^2
      if (P2 == P0) {
        P2 = coprime_near(np, 0.2360679774997897L);            // 1/phi^3
        if (P2 == P0) { P2 = 1u; }
      }

      // sigma_s = sqrt(kT_s / m_s), matching the injector's convention that
      // setup.temperatures is divided by the species mass (archetypes/utils.h:59)
      const auto m1   = domain.species[0].mass();
      const auto m2   = domain.species[1].mass();
      // The stock injector forms theta_s = temperatures[s]/mass_s and hands it
      // to JuttnerSinge; for theta << 1 (here 1e-3 and 5.4e-7) that is a
      // Gaussian of width sqrt(theta). We load that Gaussian directly, which is
      // also what SHARP does -- non-relativistic, and 1V when one_v is set.
      const bool shear = (loading == "quiet");

      // The stratified midpoint set {probit((j+1/2)/ppc)} has variance BELOW
      // unity (-1.97% at ppc=64, -0.64% at 200), so "quiet_lattice" measures the
      // realised variance and divides it out to land on the exact temperature.
      //
      // "quiet" must NOT do this. Its velocities are sampled, so the realised
      // variance carries a physical 1/sqrt(N) fluctuation -- the same sampling
      // whose density counterpart is the shot noise we are restoring. Dividing
      // it out would re-quiet the very channel this loader exists to open.
      real_t vnorm { ONE };
      if (not shear) {
        real_t sumsq { ZERO };
        Kokkos::parallel_reduce("QuietNorm", ppc, QuietNorm_kernel(ppc), sumsq);
        vnorm = math::sqrt(sumsq / static_cast<real_t>(ppc));
      }

      const auto sig1 = math::sqrt(temperatures[0] / m1) / vnorm;
      const auto sig2 = math::sqrt(temperatures[1] / m2) / vnorm;

      Kokkos::parallel_for(
        "QuietInit",
        npart,
        QuietInit_kernel<S, M>(domain.species[0], domain.species[1], xi_min, ppc,
                               ncell, P0, P2, sig1, sig2, one_v, shear));
      for (auto sp : { 0u, 1u }) {
        domain.species[sp].set_npart(domain.species[sp].npart() + npart);
        domain.species[sp].set_counter(domain.species[sp].counter() + npart);
      }
    }
  };

} // namespace user

#endif
