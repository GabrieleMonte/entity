#ifndef PROBLEM_GENERATOR_H
#define PROBLEM_GENERATOR_H

#include "enums.h"
#include "global.h"

#include "arch/kokkos_aliases.h"
#include "traits/pgen.h"
#include "utils/error.h"
#include "utils/numeric.h"

#include "framework/domain/domain.h"
#include "framework/domain/metadomain.h"
#include "framework/parameters/parameters.h"

#include <string>
#include <utility>

namespace user {
  using namespace ntt;

  enum class FieldGeometry : uint8_t {
    dipole,
    monopole
  };

  template <Dimension D>
  struct InitFields {
    InitFields(real_t bsurf, real_t rstar, const std::string& field_geometry)
      : Bsurf { bsurf }
      , Rstar { rstar }
      , field_geom { field_geometry == "monopole" ? FieldGeometry::monopole
                                                  : FieldGeometry::dipole } {}

    Inline auto bx1(const coord_t<D>& x_Ph) const -> real_t {
      if (field_geom == FieldGeometry::monopole) {
        return Bsurf / SQR(x_Ph[0] / Rstar);
      } else {
        return Bsurf * math::cos(x_Ph[1]) / CUBE(x_Ph[0] / Rstar);
      }
    }

    Inline auto bx2(const coord_t<D>& x_Ph) const -> real_t {
      if (field_geom == FieldGeometry::monopole) {
        return ZERO;
      } else {
        return Bsurf * HALF * math::sin(x_Ph[1]) / CUBE(x_Ph[0] / Rstar);
      }
    }

  private:
    const real_t        Bsurf, Rstar;
    const FieldGeometry field_geom;
  };

  /*
    Rotating-conductor surface: E = -(Omega x r) x B, with Omega along the
    polar axis. `omega_eff` is the (possibly ramped) spin -- see PGen::spin().
  */
  template <Dimension D>
  struct DriveFields : public InitFields<D> {
    DriveFields(real_t             omega_eff,
                real_t             bsurf,
                real_t             rstar,
                const std::string& field_geometry)
      : InitFields<D> { bsurf, rstar, field_geometry }
      , Omega { omega_eff } {}

    using InitFields<D>::bx1;
    using InitFields<D>::bx2;

    Inline auto bx3(const coord_t<D>&) const -> real_t {
      return ZERO;
    }

    Inline auto ex1(const coord_t<D>& x_Ph) const -> real_t {
      return Omega * bx2(x_Ph) * x_Ph[0] * math::sin(x_Ph[1]);
    }

    Inline auto ex2(const coord_t<D>& x_Ph) const -> real_t {
      return -Omega * bx1(x_Ph) * x_Ph[0] * math::sin(x_Ph[1]);
    }

    Inline auto ex3(const coord_t<D>&) const -> real_t {
      return ZERO;
    }

  private:
    const real_t Omega;
  };

  /*
    Dark-photon drive.

    The dark matter is non-relativistic (v ~ 1e-3), so A' is coherent over its
    de Broglie wavelength 2*pi/(m v), which dwarfs the resonant layer. The drive
    is therefore UNIFORM IN SPACE and oscillates at omega = m_A'; its amplitude
    e0 = eps * m_A' * A' is evaluated once per step by PGen::ext_efield().

    Polarization is along z-hat (the rotation/dipole axis). That is forced, not
    chosen: this engine is 2D axisymmetric, and an x/y-polarized drive carries a
    cos(phi) that cannot be represented. A uniform z-hat field in the ORTHONORMAL
    spherical (tetrad) basis the pusher expects -- see kernels/pushers/sr.hpp:209,
    which applies transform_xyz<Idx::T, Idx::XYZ> to (ex1, ex2, ex3) -- is

        E_rhat     =  e0 cos(theta)
        E_thetahat = -e0 sin(theta)

    NOTE these are felt by the particles but are NOT evolved by Maxwell. The
    particles respond, deposit current, and that current sources the real
    outgoing EM field. That is the conversion.
  */
  struct ExtFields {
    ExtFields(real_t e0) : e0 { e0 } {}

    /*
      Generic over the coordinate-array length, deliberately. The pusher tests
      for these methods against coord_t<M::Dim> (2D here, sr.hpp:70) but calls
      them with coord_t<M::PrtlDim>, which is 3D in spherical (sr.hpp:201,
      metrics/spherical.h:40). Minkowski hides the mismatch because PrtlDim == D
      there, so the dpdm version of this struct could hardcode coord_t<D>. A
      member template satisfies detection and call site both. x_Ph[1] is theta
      either way.
    */
    template <class C>
    Inline auto ex1(const C& x_Ph) const -> real_t {
      return e0 * math::cos(x_Ph[1]);
    }

    template <class C>
    Inline auto ex2(const C& x_Ph) const -> real_t {
      return -e0 * math::sin(x_Ph[1]);
    }

  private:
    const real_t e0;
  };

  /*
    Outer-boundary sponge target: the STATIC VACUUM DIPOLE.

    InitFields alone defines only bx1/bx2, and MatchBoundaries_kernel ignores
    any component a setter does not define (kernels/fields_bcs.hpp:38). So a
    MatchFields returning InitFields relaxes ONLY B_r and B_theta and leaves E
    and B_phi completely untouched -- i.e. there is no absorbing boundary for
    outgoing EM waves at all, at any value of grid.boundaries.match.ds.

    Naming bx3/ex1/ex2/ex3 (all zero in the static vacuum dipole) is what turns
    the layer into a real sponge. Inside it the wind's genuine B_phi is damped
    too, which is the usual price of a sponge and why nothing in that layer
    should be trusted or plotted.
  */
  template <Dimension D>
  struct MatchVacuum : public InitFields<D> {
    MatchVacuum(real_t bsurf, real_t rstar, const std::string& g)
      : InitFields<D> { bsurf, rstar, g } {}

    using InitFields<D>::bx1;
    using InitFields<D>::bx2;

    Inline auto bx3(const coord_t<D>&) const -> real_t { return ZERO; }
    Inline auto ex1(const coord_t<D>&) const -> real_t { return ZERO; }
    Inline auto ex2(const coord_t<D>&) const -> real_t { return ZERO; }
    Inline auto ex3(const coord_t<D>&) const -> real_t { return ZERO; }
  };

  template <SimEngine::type S, class M>
  struct PGen {
    static constexpr auto D { M::Dim };
    // compatibility traits for the problem generator
    static constexpr auto engines {
      ::traits::pgen::compatible_with<SimEngine::SRPIC> {}
    };
    static constexpr auto metrics {
      ::traits::pgen::compatible_with<Metric::Spherical, Metric::QSpherical> {}
    };
    static constexpr auto dimensions { ::traits::pgen::compatible_with<Dim::_2D> {} };

    const real_t      Bsurf, Rstar, Omega, spinup;
    const real_t      A0, omega_d, drive_ramp;
    const std::string field_geom;
    InitFields<D>     init_flds;

    PGen(const SimulationParams& p, const Metadomain<S, M>& m)
      : Bsurf { p.template get<real_t>("setup.Bsurf", ONE) }
      , Rstar { m.mesh().extent(in::x1).first }
      , Omega { static_cast<real_t>(constant::TWO_PI) /
                p.template get<real_t>("setup.period", ONE) }
      , spinup { p.template get<real_t>("setup.spinup", ZERO) }
      , A0 { p.template get<real_t>("setup.A0", ZERO) }
      , omega_d { p.template get<real_t>("setup.omega_d", ZERO) }
      , drive_ramp { p.template get<real_t>("setup.drive_ramp", ZERO) }
      , field_geom { p.template get<std::string>("setup.field_geometry", "dipole") }
      , init_flds { Bsurf, Rstar, field_geom } {
      raise::ErrorIf(A0 != ZERO and omega_d <= ZERO,
                     "setup.omega_d must be > 0 when the drive is on (A0 != 0)",
                     HERE);
    }

    /*
      Spin, smoothly ramped from 0 to Omega over `setup.spinup`. Switching the
      surface E on instantaneously against a vacuum dipole launches a spurious
      EM pulse; sin^2 ramps it in with zero initial slope. spinup = 0 restores
      the instantaneous switch-on.
    */
    auto spin(simtime_t time) const -> real_t {
      const auto t = static_cast<real_t>(time);
      if (spinup <= ZERO or t >= spinup) {
        return Omega;
      }
      return Omega * SQR(math::sin(static_cast<real_t>(constant::HALF_PI) * t / spinup));
    }

    /*
      A0*cos(omega_d*t) starts at FULL amplitude at t=0, which kicks every
      particle impulsively and rings the whole box at its local omega_p --
      broadband noise that buries a narrow resonance. `setup.drive_ramp` fades
      the amplitude in over that many time units (sin^2, zero initial slope).
      drive_ramp = 0 keeps the original abrupt switch-on.
    */
    auto ext_efield(simtime_t time) const -> real_t {
      const auto t   = static_cast<real_t>(time);
      auto       amp = A0;
      if (drive_ramp > ZERO and t < drive_ramp) {
        amp *= SQR(
          math::sin(static_cast<real_t>(constant::HALF_PI) * t / drive_ramp));
      }
      return amp * math::cos(omega_d * t);
    }

    /*
      @note the drive is uniform, so cos() is evaluated once per species per
            step here rather than once per particle in the kernel
      @note applied to both species; the pusher scales by (q/m) per species
    */
    auto ExternalFields(simtime_t time, spidx_t, const Domain<S, M>&) const
      -> std::pair<bool, ExtFields> {
      return {
        A0 != ZERO,
        ExtFields { ext_efield(time) }
      };
    }

    // records the instantaneous drive alongside the response
    auto CustomStat(const std::string&  name,
                    timestep_t,
                    simtime_t           time,
                    const Domain<S, M>&) const -> real_t {
      if (name == "e_ext") {
        return ext_efield(time);
      } else if (name == "omega") {
        return spin(time);
      } else {
        raise::Error("Unrecognized custom stat: " + name, HERE);
        return ZERO;
      }
    }

    auto AtmFields(simtime_t time) const -> DriveFields<D> {
      return DriveFields<D> { spin(time), Bsurf, Rstar, field_geom };
    }

    auto MatchFields(simtime_t) const -> MatchVacuum<D> {
      return MatchVacuum<D> { Bsurf, Rstar, field_geom };
    }
  };

} // namespace user

#endif
