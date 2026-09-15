# Reproducing arXiv:2510.13956 in Entity

Working plan for simulating the dark-photon-dark-matter (DPDM) resonant-conversion
saturation of Hook, Huang & Shalaby (*"No cosmological constraints on dark photon
dark matter from resonant conversion: Impact of nonlinear plasma dynamics"*,
arXiv:2510.13956v1) with **Entity** instead of their **SHARP** code.

Companion to `PIC_NOTES.md`. That file is the ladder from zero plasma background to
running Entity; this file is a single-paper reproduction plan and assumes Runs 1–3
of that ladder are understood (two-stream, Weibel, reconnection).

Paper artefacts for comparison: figures 2, 3, 6, 7, 8 and the animations at
`mohamadshalaby.github.io/dpdm.html` (their ref. [59]).

> **Status (2026-09-15).** **Both resonance cases of the paper are reproduced.**
> The quiet start is implemented, validated and in production use (§3b). The paper
> source (`arXiv-2510.13956v1.tar.gz`) has been read and every simulation parameter
> is now taken from it rather than inferred (§0e). Landau–Zener is under way (§8).
>
> Everything in §0 is *measured*, not estimated. Two earlier headline claims in this
> file were **wrong and are now corrected**:
>
> - "the shot-noise floor is ~650x below Eq. B4" — a units error. `Nₚ` in Eq. B4 is
>   the **total** macroparticle count, not per-cell. Correctly read, Eq. B4 agrees
>   with measurement to a flat 3.1x (§0c). **Use it to size `Np`.**
> - "the quiet start is very likely unnecessary" — false on both counts. It *is*
>   necessary (the paper uses one), and it does **not** remove the particle-count
>   requirement: the ordered state relaxes to the thermal shot-noise floor by
>   `ωₚt ≈ 5000` (§0f). It fixes `t = 0` only.
>
> Cost is also far above the old "~20 node-hours": the slow resonance case alone took
> **~51 GPU-hours** because per-step cost roughly doubles once the plasma saturates
> (§0g). That effect is the single most expensive lesson in this file.

---

## 0. Measured calibration (Lonestar6, `esirkepov` + `shape_order=5` + `double`)

All runs: `Nx = 1000`, `L = 40 c/ωₚ`, `CFL = 0.5` (`dt = 0.02`), `current_filters = 0`,
`skindepth0 = larmor0 = 1`, `densities = [2.0]`, `temperatures = [1e-3, 1e-3]`.
`Np` below is `ppc0`, i.e. macroparticles per cell **per species**.

### Unit conventions, pinned empirically

| Quantity | Result |
|---|---|
| `nₑ` | `T00_1` rest-mass part `= 1.0000` ⇒ `nₑ = n0` exactly ⇒ **`ωₚₑ = 1`** (the §4 √2 trap is handled correctly by `densities = [2.0]`) |
| thermal energy | `T00_1 − 1 = T00_2 − 1836 = 1.5e-3 = (3/2) kT` per species (3V Maxwellian), as expected |
| **field energy** | `ε_E / (n₀m₀c²) = ½ (E1^2 + E2^2 + E3^2)` with the stats columns **as reported** — fixed by the conservation fit `ΔT00 = −½Δ⟨E²⟩`, `r = −0.996` |
| `ωₚₑ` | detuning scan (below) puts the resonance at `ω = 1.000`, symmetric, to better than 2% |

### Energy conservation — the §9 go/no-go: **pass**

Undriven, `Np = 500`, `ωₚt = 400` (`2×10⁴` steps): after the noise plateau sets in
(`ωₚt ≈ 2`) the electron energy drifts at `−4.5e-10` per `ωₚ⁻¹`, and the loss appears
as an ion gain — noise-mediated equilibration, **not grid heating**. Extrapolated to
`ωₚt = 8×10⁴` that is 2.4% of the thermal energy, and it scales as `1/Np`, so 0.6% at
`Np = 2000`. Yee + Esirkepov survives the `4×10⁶`-step run.

### Shot-noise floor — Fig. 5, and Eq. B4 is **right** (earlier note was wrong)

```
ε_noise = ½⟨E₁²⟩ = 2.06e-4 / ppc0          [n_e m_e c²]   (ppc0 per cell per species)
```

Verified over `ppc0 ∈ {500, 2000, 8000}`: `ε_noise · ppc0` constant to 3%. Reached by
`ωₚt ≈ 2` (the co-located pair start gives `ε_noise = 8e-11` at `t = 0`, i.e. zero)
and flat thereafter.

> ⚠ **Correction (2026-09-14, from the paper source).** An earlier version of this
> section claimed Eq. B4 was "650x off" and "a bare-Poisson expression that should
> not be used to size `Np`". **That was a units error and is wrong.** Eq. B4 is
> `ε_noise = L²/(12 Nₚ Nₓ)` where the paper states explicitly that **`Nₚ` is the
> TOTAL number of macroparticles**, not the per-cell count. The 650x came from
> substituting `Nₚ = ppc0`. Read correctly (`Nₚ = ppc0 · Nₓ · 2` for two species):

| `ppc0` | Eq. B4 | measured | ratio |
|---|---|---|---|
| 200 | 3.33e-7 | 1.03e-6 | 3.1x |
| 2000 | 3.33e-8 | 1.03e-7 | 3.1x |
| 200000 | 3.33e-10 | 1.03e-9 | 3.1x |

A consistent **3.1x**, flat across three decades — i.e. Eq. B4 has the right scaling
and the right magnitude to within a small constant (plausibly a species-counting or
`⟨E²⟩`-vs-potential-energy convention). **Use it to size `Np`.** The paper's own
`N_pc` choices follow from it and should be matched rather than argued down.

The screened-form-factor calculation that "explained" the 650x was explaining an
artefact of the misreading.

### Drive normalization — validated end-to-end

`v_q^D/vₜₕᵉ = 0.03` (`A0 = 9.4868e-4`), `ω = 1`, `Np = 2000`, least squares over
`10 < ωₚt < 40`:

```
ε_E + ΔKE_e  =  C · A0² (ωₚ t)²        C = 0.12432    vs Eq. A25's 1/8 = 0.125
```

0.5%. **`A0_code = A0_paper = v_q^D/c`** — §4's recommended choice is correct as
written. Two riders:

- It is the **sum** that obeys the `t²` law. The field alone is half of it on
  cycle-average and swings between 3% and 85% of the total at `2ωₚ`. Compare Fig. 2's
  linear-theory line against `ε_E + ΔKE_e`, never `½⟨Eₓ²⟩` alone.
- Cycle-averaged equipartition `⟨ε_E⟩/⟨ΔKE_e⟩ = 1.12`.

Ion response: `ΔKE_i/ΔKE_e = 5.70e-4` vs `1/1836 = 5.45e-4` — the drive reaches both
species with the right `q/m` scaling.

### `ωₚ` calibration — detuning scan

`v_q^D/vₜₕᵉ = 10⁻³` (`A0 = 3.162e-5`), `Np = 2000`, `ωₚt = 400`, five frequencies.
`ε_tot(400)` normalized to the `ω = 1` case, against the linear-theory `sinc²` envelope
`|sin(δt/2)/(δt/2)|²`:

| `ω` | 0.96 | 0.98 | **1.00** | 1.02 | 1.04 |
|---|---|---|---|---|---|
| measured | 0.016 | 0.022 | **1.000** | 0.048 | 0.013 |
| `sinc²` | 0.015 | 0.036 | 1.000 | 0.036 | 0.015 |

`ω = 1` grows as `t²` across a factor 100 in energy; every detuned case stays bounded
and beats. Symmetric about 1.00. The on-resonance run holds `C ≈ 0.119` from
`ωₚt = 34` to `370`.

### Landau–Zener swept phase

`e_ext(t)` from the `CustomStat` matches `A0 cos(ω₀t + ½ ω̇ t²)` to **2.8e-14** (double
precision, exact) once the one-`dt` diagnostic offset below is applied. The naive
`ω(t)·t` form would deviate by up to `0.86 A0` — §6.6 was right to warn.

### Noise budget, rebuilt on the measured floor

```
t̄_noise  = 0.0406 / (A0 √Np)             (ωₚt at which ε_driven = ε_noise)
Np        = (0.0406 / (A0 t̄))²            (Np needed for a target t̄)
t_sat     = 0.0632 / A0                   (ωₚt at which ε_driven = ½ nₑ vₜₕᵉ²)
```

The useful consequence:

```
t̄_noise / t_sat = 0.642 / √Np            — independent of A0
```

At `Np = 2000` the driven signal dominates the noise over the last **98.6% of the
growth phase, at every amplitude**. The paper's `t̄ ~ 2.4` is an *absolute* time
criterion and is needlessly strict at small `A0`, where saturation takes
proportionally longer anyway. **Recommendation: size `Np` by a convergence study
(§8 run 3), not by a formula.**

For the LZ family the criterion is genuinely different, because before resonance the
drive is a *bounded* forced oscillation of amplitude `A0/(1−ω²/ωₚ²) = 2.78 A0` at
`ω = 0.8 ωₚ`, not a growing one. Requiring noise `50×` below that forced energy
(`ε_forced = 1.93 A0²`) gives `Np ≥ 5.34e-3 / A0²`:

| `v_q^D/vₜₕᵉ` | 1.0 | 0.1 | 0.03 |
|---|---|---|---|
| `Np` required | 5 | 535 | 5930 |

⚠ **Superseded (2026-09-14).** This table was computed from the measured floor
*after* wrongly discarding Eq. B4 (see the correction above). The paper sizes the LZ
runs from its own stated criterion — `t̄_noise ≤ 0.02`, i.e. driven energy 50x shot
noise — which with `Nₚ` read correctly as the total macroparticle count gives

| `v_q^D/vₜₕᵉ` | 1.0 | 0.1 |
|---|---|---|
| `N_pc` required | `2.7e3` | `2.7e5` |

These are the numbers to use. The `10³` figure below is not reliable.

### Cost model

```
2.42 ns per particle-step per node   (1 LS6 CPU node, 112 threads, Np ≥ 8000)
node-hours = 6.7e-8 × Np × (ωₚ t)    for Nx = 1000
```

OpenMP scaling 14 → 112 threads: 5.83x (73%). Deposit uses `Kokkos::ScatterView`, so
the 1000-cell grid causes no atomic contention on CPU.

---

### 0e. Paper parameters, taken from the source

`arXiv-2510.13956v1.tar.gz` contains `main.tex`. These are quoted, not inferred:

| quantity | value | where |
|---|---|---|
| grid | `Nx = 1000`, `L = 40 c/ωₚ`, `Δx = 0.04 c/ωₚ` | "In all simulations, we use 1000 computational cells" |
| temperatures | `kᵦTₑ = kᵦTᵢ = 10⁻³ mₑc²`, `λ_D ≈ 0.033 c/ωₚ` | Eq. dist |
| **`vₜₕ` definition** | `vₜₕ² ≡ ⟨(vₓ − ⟨vₓ⟩)²⟩ = kᵦT/m` — **1D, x-component only** | after Eq. dist |
| thermal energy density | `5×10⁻⁴` in `nₑmₑc²` | §"validate the shot noise" |
| resonance `N_pc` | `2×10²` and `2×10⁵` | Fig. A00 |
| LZ sweep | `ω(t) = 0.8ωₚ + 0.1ωₚ t/5000` ⇒ `dω/dt = 2.0e-5` | app. sim-params |
| LZ amplitudes | `v_q^D/vₜₕᵉ = 1.0` (top) and `0.1` (bottom) | Fig. LZcase caption |
| LZ `N_pc` criterion | `t̄_noise ≤ 0.02`, i.e. driven energy 50x shot noise | app. sim-params |
| shot noise | `ε_noise = L²/(12 Nₚ Nₓ)`, **`Nₚ` = TOTAL macroparticles** | Eq. PE_noise |

The `vₜₕ` definition independently vindicates the **1V** choice (`one_v = true`,
`ux2 = ux3 = 0`) and the `ux1`-variance temperature in §6.5. The quoted `5×10⁻⁴`
matches our `eps_thermal_1V` reference exactly.

LZ crossing is at `ωₚt = 10⁴` (`0.8 + 2e-5·t = 1`). The paper's "nonlinearity kills
the resonant transfer around `ωₚt = 5000`" refers to `ω = 0.9ωₚ`, i.e. **before** the
crossing — off-resonant excitation, not the crossing itself.

Solving `t̄ = L/(A0√(3NₚNₓ/2)) ≤ 0.02` for the LZ runs:

| `v_q^D/vₜₕᵉ` | `A0` | `N_pc` |
|---|---|---|
| 1.0 | 3.16228e-2 | 2.7e3 |
| 0.1 | 3.16228e-3 | 2.7e5 |

### 0f. The quiet start relaxes — it fixes `t = 0` only

Undriven, `ppc0 = 2700`, `A0 = 0`, `ωₚt = 2×10⁴` (10⁶ steps), quiet loading:

```
ε_E:  6.7e-23  (t = 0)   →   4.8e-8  plateau, reached by ωₚt ≈ 5000
random-start floor at the same ppc0                = 7.6e-8
```

So the ordered lattice starts at round-off and **relaxes to ~62% of the thermal
shot-noise floor**. This is the paper's own statement ("it takes approximately
`tωₚ ~ 3–4` for the driven energy to be converted into de-coherent high-`k` noise")
seen from the undriven side. **Consequence: a quiet start does not let you reduce
`Np`.** Size `Np` from Eq. B4 / the paper's criterion; the quiet start buys a clean
`t = 0` and nothing more.

Same run, numerical heating over 10⁶ steps:

```
kTₑ drift = −0.0162%      kTᵢ drift = +0.0084%
```

Negligible. An ordered lattice does **not** drift — this was the gate for LZ, where a
slow `Tₑ` creep would detune `ωₚ` and masquerade as physics.

### 0g. Per-step cost doubles after saturation — the expensive lesson

Measured on one A100 at `ppc0 = 2×10⁵` (4×10⁸ particles):

```
265 ms/step   laminar, pre-saturation
516 ms/step   after the plasma saturates and goes turbulent
```

A 2000-step benchmark (`ωₚt = 40`) showed a **dead-flat 265 ms** and was taken as
proof of stability. It was nothing of the kind: the run saturates at `ωₚt ≈ 2000`,
so the benchmark never left the laminar phase. The slow resonance case consequently
overran its 48 h wall at 88% and needed a checkpoint restart.

**Rule: benchmark across the physics transition, or assume 2x.** The CPU showed the
same effect (+48%, plateauing by step ~300 in the driven case) and it was wrongly
dismissed as a host-memory artefact that would not apply on GPU.

### 0h. Machine performance summary (measured)

| configuration | ms/step @ `ppc0 = 2e5` | note |
|---|---|---|
| 1 CPU node, 128 cores | ~950 (steady state) | rises +48% from 640 |
| 8 CPU nodes, MPI | 135 | 4.94x speedup, 69% efficiency |
| **1× A100** | **265 → 516** | 3.6x a full CPU node |
| 3× A100, 1 node | 128 | 2.07x over 1 GPU, 69% efficiency |

GPU deposition atomics on a 1000-cell grid cost only **1.42x** going from
`ppc0 = 2e4` to `2e5` — the feared `ScatterAtomic` collapse did not happen.

On a busy machine **queue wait dominates**: 6- and 12-GPU benchmark jobs sat 30 h
unscheduled and were cancelled, while 1-node jobs ran within minutes. Prefer fewer
nodes than the scaling curve suggests.

---

## 1. The paper as a PIC problem

Strip the cosmology and it is a very small kinetic run.

| Ingredient | Paper's choice |
|---|---|
| Geometry | 1D-1V **electrostatic**, periodic box |
| Grid | `Nx = 1000` cells, `L = 40 c/ωₚ`, `Δx = 0.04 c/ωₚ` |
| Species | electrons + ions, **real mass ratio `mᵢ/mₑ = 1836`**, both evolved self-consistently |
| Init spatial | uniform density, electrons and ions **co-located**, particles **equally spaced** (quiet start) so that `ε_noise = 0` at `t = 0` |
| Init velocity | isotropic Maxwellian, `kᵦTₑ = kᵦTᵢ = 10⁻³ mₑc²` ⇒ `vₜₕᵉ = 3.1×10⁻² c`, `vₜₕⁱ = 7.4×10⁻⁴ c` |
| Debye length | `λ_D ≈ 0.033 c/ωₚ` — deliberately `≈ Δx`, "marginally resolved", to suppress grid heating (Birdsall & Maron 1980) |
| Driver | the dark photon enters as an **external, spatially uniform (`k = 0`), time-oscillating electric field** `E_ext(t) = A₀ cos(ωt)`, with `A₀ ≡ εeE′/(mₑωₚc) = v_q^D/c`. It acts on the particle Lorentz force only; it is **not** a source term in Maxwell's equations. |
| Pusher | Vay (2008) |
| Deposition / interpolation | **5th-order splines**, charge + current; total momentum conserved exactly; energy "well controlled" |
| Timestep | `dt ≈ 0.02 /ωₚ` (`Δx = 0.04`, CFL `= 0.5`) |
| Longest run | `ωₚt ≈ 8×10⁴` (≈ `4×10⁶` steps) |

### Two run families

- **Resonance:** `ω = ωₚ` fixed. Cleanest illustration; more efficient than the
  realistic case.
- **Landau–Zener (LZ):** `ω(t)` ramped slowly through `ωₚ`
  (`ω(t) = 0.8 ωₚ + 0.1 ωₚ · t/5000`, i.e. a slow linear chirp — a numerically
  tractable proxy for the real `ωₚ(t)` drift, since `ωₚ/H ~ 10¹⁰` periods is
  unsimulable).

### Scanned parameter

`v_q^D / vₜₕᵉ` from `10⁻³` to `~1`, i.e. `A₀` from `~3×10⁻⁵` to `~3×10⁻²`.
Two qualitative regimes:

- **slow growth** (`v_q^D/vₜₕᵉ ≲ (mₑ/mₚ)^{1/2}/2 ≈ 0.01`): resonance slows and
  saturates *before* the Langmuir energy reaches the electron thermal energy;
- **fast growth** (`v_q^D/vₜₕᵉ ≳ 0.01`): electric-field energy overshoots the
  thermal energy, then halts after `~1/ωₚᵢ` when the ions respond; electrons end
  `O(10²)` times hotter, `Tₑ/Tᵢ → O(30)`, steady state.

### The physics being demonstrated

The `k = 0` Langmuir wave grows linearly (`E ∝ t`, energy `∝ t²`, App. A Eq. A25)
until the ponderomotive force

```
F_p = −∇Φ_p ,   Φ_p(x) = e² Ê²(x) / (4 mₑ ω²)
```

pumps `k ≠ 0` Langmuir waves (modulational instability) and ion-acoustic waves
(electrostatic-decay / ion-acoustic instability). Those make `ntot(x)` and hence
`ωₚ(x)` spatially non-uniform, which detunes the resonance and shuts conversion
off once the deposited energy `≈` the electron thermal energy. Ion-acoustic waves
then damp only on the slow electron–ion collision timescale (`τei`), which is
collisionless-frozen in the PIC, so the suppressed state persists.

Net cosmology result (not simulated, downstream of the PIC): DPDM `ε` limits from
resonant conversion weaken by `3000` to `10⁷` over ten decades in `mₐ′`.

---

## 2. What Entity already provides

Most of the machinery exists. File references are to this repo.

| Need | Entity mechanism |
|---|---|
| 1D SRPIC / Minkowski | registered specialization, `src/framework/specialization_registry.h:41`; identical class of setup to `pgens/streaming/twostream.toml` |
| longitudinal `Eₓ` from deposited `Jₓ`, `B` stays ~0 | two-stream "free sanity check" in `PIC_NOTES.md §2` |
| **external, time-varying, `k=0` E field on particles only** | `PGen::ExternalFields(time, sp, domain)` → `{apply?, ExtFields}`; added to interpolated grid E **for the force only** in `src/kernels/pushers/sr.hpp:197-213`; never enters deposit or Ampère. Demonstrated by `examples/external_fields/`. |
| species-dependent drive | `ExternalFields` receives `sp`; pusher multiplies by `q/m` per species (`normalized_dt_half`, `src/kernels/pushers/sr.hpp:111`) |
| arbitrary mass ratio | `mass = 1836.0` on the ion species (TOML) |
| Maxwellian init at fixed `T`, co-located pairs | `arch::InjectUniformMaxwellians`, `src/archetypes/utils.h:47`; pair species injected at the **same** cell position, `src/kernels/injectors.hpp:251-286` ⇒ `ρ(x, t=0) = 0` exactly |
| 5th-order splines | `-D deposit=esirkepov -D shape_order=5` (CMake allows 1–11, `CMakeLists.txt:94`; `S5` in `src/kernels/particle_moments.hpp:124`) |
| Vay pusher | `velocityEMPush_Vay`, `src/kernels/pushers/sr.hpp:370` |
| per-species kinetic-energy time series | `[output.stats]` `T00_1`, `T00_2` (bare `T00` gives the **summed** value over species, `= 1837.003` here). `T00` moment `= Σ wₚ mₚ γₚ` (`src/kernels/particle_moments.hpp:214`); rest-mass part constant ⇒ `ΔT00_s = ΔKE_s` |
| EM field energy, Poynting, work | `[output.stats]` `E^2`, `B^2`, `ExB`, **`J.E`** (`src/global/enums.h:286`). ⚠ the enum lookup is lowercase `"j.e"` — `"JdotE"` **aborts** at `src/output/stats.cpp:42`. `E^2` splits into `E1^2`/`E2^2`/`E3^2` columns |
| real-space panels | `[output.fields]` `N_1`, `N_2`, `Rho`, `E`, `J`, `T00_1`, `T00_2` (`FldsID`, `src/global/enums.h:248`) |
| checkpoint/resume for `4×10⁶`-step runs | on by default; `walltime` guard under `[checkpoint]` (`PIC_NOTES.md §6`) |

### Unit system recap (`src/framework/parameters/parameters.cpp:47-78`)

- `skindepth0` → `ωₚ0 = 1/skindepth0` for the **reference** density `n0`.
- `larmor0` → `B0 = ωB0 = 1/larmor0`. No background B here, so `larmor0` is a pure
  normalization knob.
- `sigma0 = (skindepth0/larmor0)²` — irrelevant (no B), leave it wherever.
- `n0 = ppc0 / V0`, `q0 = V0 / (ppc0 · skindepth0²)`.
- **Field vs particle stats normalization** (measured, §0): the `E^2` columns and
  `T00_*` are in the same units, with `ε_E = ½(E1^2+E2^2+E3^2)` — no `sigma0` factor
  needed at `skindepth0 = larmor0 = 1`.

---

## 3. What had to be built — **done**, see `pgens/dpdm/`

### 3a. A dedicated pgen — `pgens/dpdm/` ✅ written

`pgens/dpdm/pgen.hpp` (163 lines) plus `resonance.toml` and `landau_zener.toml`.
Auto-discovered by the CMake glob (`cmake/config.cmake:33`) — no registry change
needed, 1D SRPIC Minkowski is already in `NTT_FOREACH_SPECIALIZATION`.

Two implementation notes beyond the sketch below:

- the drive is uniform in space, so `cos()` is evaluated **once per species per step
  on the host** inside `ExternalFields` and only the resulting scalar is handed to
  `ExtFields`, rather than once per particle in the kernel;
- `ExtFields` deliberately defines **only** `ex1`, so the `HasExtEx2`/`HasExtBfield`
  branches of `src/kernels/pushers/sr.hpp` compile out entirely.

- **`InitPrtls`** — two species (e⁻: `mass 1`, `charge -1`; ion: `mass 1836`,
  `charge +1`), one uniform Maxwellian pair at `T = 10⁻³`. Reuse
  `arch::InjectUniformMaxwellians`.
- **`ExternalFields(time, sp, domain)`** — return `{ true, ExtFields{...} }` for
  **every** `sp` (drive both; ions receive `1/1836` of the kick automatically).
  `ExtFields::ex1()` returns:
  - resonance: `A0_code * cos(omega * time)`;
  - LZ: `A0_code * cos(phi(time))` with the **analytically integrated** swept
    phase `phi(t) = ω₀ t + ½ (dω/dt) t²`. Do **not** use `ω(t)·t` — it accumulates
    a spurious chirp.
- **compat traits** — `SimEngine::SRPIC`, `Metric::Minkowski`, `Dim::_1D`
  (optionally also `_2D`/`_3D` to verify the paper's claim that 1+1D captures the
  physics — see §8).
- **`CustomStat`** — implemented, exposes `e_ext` to the stats CSV via
  `[output.stats] custom = ["e_ext"]`, so the drive and the response sit in the same
  file. ⚠ it receives `time` **one `dt` behind** the CSV's own `time` column
  (`engine.hpp:311-323` passes `time - dt` to the writer); shift by one step before
  comparing. Diagnostic only — the pusher's centering is correct (§6.7).
- optional **`CustomFieldOutput`** — dump `Φ_p ∝ Eₓ²` for the animations. Not
  implemented; trivially recovered post-hoc from the `E` field dumps.

Config knobs to read from `[setup]`: `A0`, `omega`, `drive` (`"resonance"` /
`"landau_zener"`), `omega0`, `domega_dt`.

### 3b. Quiet start — **implemented and validated**, `setup.loading = "quiet"`

Entity's stock injector places particles at **random** sub-cell positions
(`src/kernels/injectors.hpp:202`). `pgens/dpdm/pgen.hpp` now adds a deterministic
path, selected by `setup.loading = "quiet" | "random"` (random stays the default so
earlier runs remain reproducible).

**The construction.** Index by velocity **class**, not by sub-cell slot:

```
c    = p / ncell                      velocity class, exactly ncell particles each
i    = p % ncell                      cell
phi  = ((c*P0) mod ppc + 1/2)/ppc     sub-cell offset — a function of c ALONE
x    = xi_min + i + phi
u_e  = sig1 * probit((c + 1/2)/ppc)
u_i  = sig2 * probit(((c*P2) mod ppc + 1/2)/ppc)
```

Four properties, each doing a specific job:

- **`phi` depends only on `c`** ⇒ every velocity class is a perfect `Δx`-spaced
  lattice spanning the box. Under streaming a class merely translates, and such a
  lattice deposits exactly uniform density at *any* translation. **This is what makes
  the order survive.** Keying the velocity to the sub-cell slot instead (with a
  per-cell rotation) fails precisely here — a given velocity then sits at a different
  offset in each cell, classes shear apart within one step, and the measured result
  was `δn/n = 6.3e-3` after a *single* step with a floor only 2x below random.
- **`P0`, `P2` coprime to `ppc`** ⇒ the index maps are bijections, so each cell
  carries the complete stratified quantile set exactly once (zero net drift, exact
  temperature per cell), while consecutive slots get *non*-consecutive quantiles so no
  cell holds a systematic `v(x)` ramp. Both strides must be coprime to `ppc` itself.
- **ions reuse the electron positions** ⇒ charge density vanishes bit-for-bit at
  `t = 0`; their velocities come through a second stride so the species are not
  phase-space correlated.
- **nothing references the cell index or domain offset** ⇒ the loading is identical
  however MPI decomposes the domain (verified, see below).

**Velocity normalization.** The stratified set `{probit((j+½)/ppc)}` has mean 0 by
symmetry but variance *below* 1 — midpoint quantiles under-sample the Gaussian tails
(−1.97% at `ppc = 64`, −0.64% at `ppc = 200`). A device-side `parallel_reduce`
measures the realized variance and divides it out, so the loaded temperature is
exact. Without it the plasma runs systematically cold at production `ppc0`.

`probit` is Acklam's rational approximation plus one Halley step against `erfc`;
verified accurate to ~1e-16 over `q ∈ [1e-8, 1−1e-8]`.

**Validation gates** (`pgens/dpdm/validate_quiet.sh` + `check_quiet.py`, `ppc0 = 64`):

| gate | result |
|---|---|
| lattice uniformity | `δn/n = 0.000e+00` — *exactly* zero |
| charge neutrality | `max|Charge| = 0.000e+00` |
| undriven noise floor | `1.31e-17` vs `3.03e-6` random — **2.3×10¹¹ lower** |
| no coherent seeded mode | loudest mode 10.6x median, power ~1e-38 |
| drive normalization | `C = 0.12531` vs `1/8` |
| under 8-domain MPI | identical to 1 domain to `1.5e-12` |

Two traps worth recording. `Rho` in Entity is **mass** density (`1 + 1836 = 1837`);
the charge density is the separate `Charge` field. And `Ex(t=0) = 0` proves nothing —
Entity initializes `E` to zero identically and *both* loadings co-locate the pair. The
discriminator is the **total** density `nₑ + nᵢ`, which is the channel the quiet start
exists to kill.

**What it does not buy.** See §0f: the ordered state relaxes to the thermal shot-noise
floor by `ωₚt ≈ 5000`. Size `Np` from Eq. B4 regardless.

### 3c. Configure-time flags — four builds, and the two traps

```sh
M="module load gcc/13.2.0"                       # the binary needs GCC 13's libstdc++
IGN=-DCMAKE_IGNORE_PATH=/work/09218/gab97/ls6/entity/cmake
COMMON="-D pgen=dpdm -D deposit=esirkepov -D shape_order=5 -D precision=double"
KOK="-D Kokkos_ENABLE_OPENMP=ON -D Kokkos_ENABLE_SERIAL=OFF"

# 1. CPU serial (OpenMP)              — validation, fast resonance case
cmake -B build-dpdm       $COMMON $KOK -D Kokkos_ARCH_ZEN3=ON \
      -D adios2_DIR=$PWD/build/_deps/adios2-build
# 2. CPU + MPI                        — multi-node fallback
cmake -B build-dpdm-mpi   $COMMON $KOK -D Kokkos_ARCH_ZEN3=ON -D mpi=ON $IGN
# 3. CUDA, single GPU                 — the production path
cmake -B build-dpdm-cuda  $COMMON -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ARCH_AMPERE80=ON $IGN
# 4. CUDA + MPI, multi-GPU            — LZ 0.1, long-term stability
cmake -B build-dpdm-cuda-mpi $COMMON $KOK -D mpi=ON -D gpu_aware_mpi=OFF \
      -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ARCH_AMPERE80=ON $IGN
```

- ⚠ **`CMAKE_IGNORE_PATH` is required** for any build that must *fetch* ADIOS2.
  `find_package(adios2 QUIET)` otherwise matches the repo's own
  `cmake/adios2Config.cmake` — which only sets cache variables and exports no
  version or targets — declares `adios2_FOUND`, and skips the fetch. Symptom:
  `ADIOS2: v` with the path pointing at `entity/cmake`, then
  `Target "entity.xc" links to adios2::cxx_mpi but the target was not found`.
  Upstream bug: *the shim shadows the package it exists to configure.*
- ⚠ **But `CMAKE_IGNORE_PATH` also blocks the Kokkos settings**, so
  `Kokkos_ENABLE_OPENMP=ON -D Kokkos_ENABLE_SERIAL=OFF` must be passed **explicitly**.
  Freshly-fetched Kokkos defaults to **Serial only**. A Serial build is
  indistinguishable from the outside and runs **~112x slower** (`NLWP = 2`, 99.5%
  CPU). This cost a 40-minute walltime and two wrong diagnoses (SLURM `-c`, then
  `ibrun` binding — both innocent). **Always verify:**
  `grep -E "^Kokkos_ENABLE_(OPENMP|SERIAL|CUDA):BOOL" build-*/CMakeCache.txt`
- **ADIOS2 2.11 renamed the C++ target** `cxx11` → `cxx` (`cxx11_mpi` → `cxx_mpi`).
  Entity hardcodes the 2.11 names (`CMakeLists.txt:161-163`), so TACC's
  `adios2/2.10.2` module *cannot* link even though it has `ADIOS2_HAVE_MPI`. Fetching
  2.11.0 (Entity's pinned tag) is the clean route.
- **`gpu_aware_mpi=OFF`** for the CUDA+MPI build: GPU-aware MPI on LS6 needs
  `mvapich2-gdr-cuda12`, a different stack from the `impi` in the gcc13 hierarchy. The
  halo here is a few ghost cells, so host staging costs nothing.
- **`esirkepov` + `shape_order=5`** — the paper's 5th-order splines; `zigzag` resets
  shape order to 1 (`CMakeLists.txt:89`).
- **`precision=double` is not optional.** `~10⁶` steps of a fast `ωₚ` oscillation
  carrying a `~10⁻⁸`-level driven signal will not survive single precision.

### 3d. Post-processing (nt2py + numpy)

- **`E(k,t)`, `n(k,t)`** — FFT the 1D field dumps. Entity's `[output.spectra]` is
  `dN/dγ`, **not** a field `k`-spectrum, so this is not built in.
- **`Tₑ(t)`, `Tᵢ(t)`** — variance of `ux1` per species from particle dumps, with
  the `k=0` bulk removed (paper's definition, their footnote 7).
- **energy-budget / rate plots** — straight from the `stats` `T00_1`, `T00_2`,
  `E^2` time series.

---

## 4. Code-units translation

Goal: pick `[scales]` and `A0_code` so numbers map cleanly onto the paper.

### Derivation of the drive amplitude

Full-step velocity kick from E in the pusher (two half-kicks of
`normalized_dt_half = ½ (q/m) ωB0 dt`, `src/kernels/pushers/sr.hpp:111`):

```
Δuₓ / dt = (q/m) · ωB0 · e_ext(t) ,   e_ext(t) = A0_code cos(Ω t)
```

For an electron (`q/m = −1`) the steady quiver solution has amplitude

```
|u_q| = ωB0 · A0_code / Ω
```

We want `|u_q| = v_q^D/c ≡ A0_paper`, driving on resonance `Ω = ωₚₑ`.

### Recommended choice

```
skindepth0 = 1.0
larmor0    = 1.0                       # ⇒ ωB0 = 1
```

Arrange the **electron** density to equal `n0` so that `ωₚₑ = 1/skindepth0 = 1`.
Mind the `PIC_NOTES.md` √2 trap: in the streaming-style pgen `densities` is
**per pair** and splits half to each species (`pgens/streaming/pgen.hpp:109-136`),
so a single `e`–`ion` pair with `densities = [1.0]` gives `nₑ = 0.5 n0` and
`ωₚₑ = 1/√2`. Use `densities = [2.0]` (⇒ `nₑ = n0`) **or** rescale
`skindepth0 = 1/√2`, then **verify**:

```python
d.fields.N_1.isel(t=0).mean()   # must be 1.0 for the recommended setup
```

With `skindepth0 = larmor0 = 1`, `nₑ = n0`, `Ω = ωₚₑ = 1`:

> **`A0_code = A0_paper = v_q^D/c = (v_q^D/vₜₕᵉ) · √(10⁻³) ≈ 0.0316 · (v_q^D/vₜₕᵉ)`**

| `v_q^D/vₜₕᵉ` | `A0_code` | regime |
|---|---|---|
| `10⁻³` | `3.2×10⁻⁵` | slow growth (linear-theory check) |
| `3×10⁻³` | `9.5×10⁻⁵` | slow growth |
| `0.03` | `9.5×10⁻⁴` | fast growth (Fig. 2 upper) |
| `0.1` | `3.2×10⁻³` | fast growth / LZ (Fig. 3 lower) |
| `1.0` | `3.2×10⁻²` | LZ (Fig. 3 upper) |

### Landau–Zener frequency

`ω(t) = ω₀ + (dω/dt) t` with `ω₀ = 0.8`, `dω/dt = 0.1/5000 = 2×10⁻⁵` (code units,
`ωₚ = 1`). Drive phase in `ex1()`:
`phi(t) = ω₀ t + ½ (dω/dt) t²`.

### Linear-theory sanity check (do this before trusting anything)

Paper App. B / Eq. A25: `ε_driven = (A0²/8)(ωₚt)² nₑ mₑ c²`. In code units
(`nₑ = n0`, `mₑ = 1`, `c = 1`), and **including the kinetic part**:

```
½⟨Eₓ²⟩ + ΔKE_e  =  (A0_code² / 8) · (ωₚ t)²
```

**Measured `C = 0.12432` vs `1/8`, i.e. 0.5%** (§0). ⚠ The earlier form of this note
had `½⟨Eₓ²⟩` alone on the left, which is wrong by a factor 2: the resonantly driven
Langmuir wave equipartitions, so the field carries only half the wave energy on
cycle-average — and instantaneously it swings between 3% and 85% of the total at
`2ωₚ`. Only the **sum** is a clean `t²`. Compare Figs. 2–3's "purple line" against
`ε_E + ΔKE_e`.

### Temperatures

`temperatures = [1e-3, 1e-3]` in units of `m0 c²`. The injector divides by species
mass (`src/archetypes/utils.h:59`) ⇒ `kᵦTₑ = kᵦTᵢ = 10⁻³ mₑc²`,
`vₜₕᵉ = √(10⁻³) c = 0.0316 c`, `vₜₕⁱ = √(10⁻³/1836) c = 7.4×10⁻⁴ c`. Matches the
paper. Verify from the `t=0` particle dump.

---

## 5. Example input — scaled-down resonance test

```toml
[simulation]
  name    = "dpdm_res"
  engine  = "srpic"
  runtime = 3000.0                     # ωₚ t

[grid]
  resolution = [1000]
  extent     = [[0.0, 40.0]]           # Δx = 0.04 = 0.04 c/ωₚ

  [grid.metric]
    metric = "minkowski"

  [grid.boundaries]
    fields    = [["PERIODIC"]]
    particles = [["PERIODIC"]]

[scales]
  skindepth0 = 1.0
  larmor0    = 1.0

[algorithms]
  current_filters = 1                  # KEEP LOW — see §6

  [algorithms.timestep]
    CFL = 0.5                          # dt ≈ 0.02 /ωₚ

[particles]
  ppc0 = 4000.0                        # scaled test; production resonance ~1e5

  [[particles.species]]
    label    = "e-"
    mass     = 1.0
    charge   = -1.0
    maxnpart = 5e6

  [[particles.species]]
    label    = "i+"
    mass     = 1836.0
    charge   = 1.0
    maxnpart = 5e6

[setup]
  temperatures = [1.0e-3, 1.0e-3]
  drive        = "resonance"           # or "landau_zener"
  A0           = 9.5e-4                # = (v_q^D/vth^e) * sqrt(1e-3)
  omega        = 1.0                   # resonance
  omega0       = 0.8                   # LZ only
  domega_dt    = 2.0e-5               # LZ only

[output]
  interval_time = 25.0

  [output.fields]
    quantities = ["E", "N_1", "N_2", "Rho", "J", "T00_1", "T00_2"]

  [output.particles]
    species = [1, 2]
    stride  = 200

  [output.stats]
    enable        = true
    interval_time = 2.0                # cheap scalars; sample the envelope
    quantities    = ["E^2", "T00_1", "T00_2", "J.E"]   # NOT "JdotE"/"T00" -- see section 2

  [output.spectra]
    enable = false

[checkpoint]
  interval_time = 500.0
  keep          = 2
```

Notes:

- `decomposition = [-1]` if you add it (single domain, no MPI) — a positive entry
  aborts in `tools::Decompose` (`PIC_NOTES.md §6`).
- `maxnpart` is a hard cap; per-species count is `ppc0 × density × 0.5 × ncells`
  (`src/archetypes/particle_injector.h:136`). Here `4000 × 2.0 × 0.5 × 1000 = 4×10⁶`
  per species — set `maxnpart` above that.
- One directory per run: you will be scanning `A0`.

---

## 6. Numerical requirements and traps

1. **Numerical heating is the #1 threat.** The long-term-stability claim runs to
   `ωₚt = 8×10⁴` (`~4×10⁶` steps). `zigzag` / 1st-order will heat the plasma past
   the physical signal long before then. Mandatory: `esirkepov` + `shape_order`
   4–5, `Δx ≲ λ_D`, `precision=double`, **and an undriven (`A0 = 0`) calibration
   run first** to measure `dTₑ/dt` and confirm it is negligible vs. the physical
   rates. This doubles as the paper's Fig. 5 / Eq. B4 shot-noise test. ✅ done; and
   repeated with quiet loading over 10⁶ steps in §0f (`kTₑ` drift −0.016%).
2. **Do not over-filter.** The physics *lives* at high `k` (modulational +
   ion-acoustic daughters up to `kλ_D ~ 1`). Entity's `current_filters` is
   binomial `[1,2,1]` smoothing (`src/kernels/digital_filter.hpp`); the 4 passes
   used in the streaming tomls will damp exactly the modes that cause saturation.
   Keep it at **0–2** and rely on the high-order shape, as the paper does.
3. ~~**`ωₚ` calibration.**~~ **CLOSED.** `skindepth0 = 1`, `densities = [2.0]` gives
   `nₑ = n0` exactly ⇒ `ωₚₑ = 1`; confirmed by a five-point detuning scan symmetric
   about `ω = 1` against the `sinc²` envelope (§0).
4. **Ion timescale.** `ωₚᵢ = ωₚₑ/√1836 ≈ 0.023 ωₚₑ` ⇒ ion-acoustic response time
   `~1/ωₚᵢ ≈ 43 ωₚₑ⁻¹`; every run must be `≫` that (they are — the fast-growth
   regime halts "after `~1/ωₚᵢ`").
5. **3V vs. 1V.** Entity carries `ux1, ux2, ux3` and its Maxwellian is isotropic;
   SHARP is genuinely 1V. In a 1D electrostatic run there are no `y,z` forces, so
   the perpendicular DOF just carry constant energy — they cancel in `ΔT00` and do
   not touch the `x` distribution that sets Landau / ion-acoustic damping. Fine,
   but: define `Tₑ` from `ux1` variance only (paper footnote 6), and do not
   compare total `T00` to the paper's 1V thermal-energy reference line without
   subtracting the inert `⊥` part.
6. **LZ phase.** Integrate the swept frequency into a phase analytically in the
   pgen (§4). Sampling `ω(t)` per step and multiplying by `t` is wrong.
7. ~~**`ExternalFields` sampling time.**~~ **CLOSED.** `pusher_ctx.time` is `tⁿ`
   (`engine.hpp:284` increments after `step_forward`), exactly the time of the grid
   `E` used in the kick. `E_ext(tⁿ)` is correctly centred; no `dt/2` correction. The
   `CustomStat` `e_ext` column receives `time - dt`, a diagnostic offset only.
8. **External field is a pure particle force, not a Maxwell source.** Correct for
   this problem (the DPDM is a `k=0` current source acting on charges) and matches
   SHARP. All feedback is through deposited `J`; confirm the 1D Ampère update
   retains `∂ₜEₓ = −Jₓ` with no transverse term (it does — two-stream confirms
   longitudinal `Eₓ` growth).
9. **Operational** (`PIC_NOTES.md §6`): one directory per run; always read back a
   `sed -i` edit before launching; `tar` before `scp`; checkpoint `walltime`
   guard for queued jobs.
10. **`ppc0` must be a TOML float.** `ppc0 = 2000` parses as integer and Entity
    aborts with `toml::value::as_floating(): bad_cast`. Write `2000.0`. Cost two
    separate job failures.
11. **Entity's default logging is two firehoses.** `<name>.log` is a VERB kernel
    trace at ~4 kB *per step* (`[diagnostics] log_level = "WARNING"`);
    `<name>.out` is a progress bar redrawn every step (`interval = 200`,
    `colored_stdout = false`). Defaults wrote 107 MB of text beside 75 kB of science
    on a 2×10⁴-step run, synchronously to the parallel filesystem.
12. **`[output.spectra]` defaults to ON.** Disable it explicitly unless wanted.
13. **Entity races on output-directory creation under MPI.** Every rank calls mkdir
    on `<simulation.name>/`; the losers abort with
    `filesystem error: cannot create directory: File exists`, exit 255, in seconds —
    which looks deceptively like a fast run. Pre-create `<rundir>/<name>/`.
14. **Checkpoint domain count is locked.** An N-domain checkpoint can only be resumed
    by N domains. Decide the GPU/rank count *before* the first submission; a
    single-GPU checkpoint cannot be resumed by a 3-GPU MPI binary.
15. **`blocking_timers = true` is the first move on any performance surprise.** It
    prints a per-kernel breakdown with ns-per-particle and answered in 9 seconds what
    two rounds of speculation got wrong.

---

## 7. Diagnostics mapping

| Paper figure | Entity source |
|---|---|
| Fig. 2 / 3 energy budget | `stats` `E^2`, `T00_1`, `T00_2`. ⚠ the `t²` line `A0²/8·(ωₚt)²` must be compared against **`ε_E + ΔKEₑ`**, not `½⟨E²⟩` alone (factor 2, §0). Thermal ref is the **1V** `5×10⁻⁴`, not `T00₁−1 = 1.5×10⁻³` (3V) |
| Fig. 6 injection rate `∂ₜE_tot/(ωₚ E_tot)` | finite-difference `T00_1 + T00_2` from `stats` |
| Fig. 7 `Tₑ/Tᵢ(t)` | `ux1` variance per species from particle dumps (bulk removed) |
| Fig. 8 real-space `nₑ, nᵢ, ρ, uₑ, uᵢ` | `fields` `N_1`, `N_2`, `Rho`; bulk `u` from `V` moment or `T0i/T00` |
| Fig. 8 bottom `E(k)`, `n(k)` | **post-hoc FFT** of `fields` `E`, `N_*` dumps (not built in) |
| Fig. 5 shot noise / Eq. B4 | undriven runs at 2–3 `ppc0`; `½⟨E²⟩` plateau vs. `L²/(12 Nₚ Nₓ)` with **`Nₚ` = TOTAL macroparticles** (= `ppc0·Nₓ·2`). Agrees to 3.1x — §0c |

Analysis skeleton:

```python
import nt2
d = nt2.Data("dpdm_res")
# energy budget
E2   = d.stats["E^2"]          # ∝ EM field energy
KEe  = d.stats["T00_1"]        # electron total energy; ΔKEe = KEe - KEe.isel(t=0)
KEi  = d.stats["T00_2"]
# linear-theory line
A0 = 9.5e-4
t  = d.stats.t
lin = A0**2 / 8 * (t)**2       # ωₚ = 1
# k-spectrum
import numpy as np
Ex  = d.fields.Ex.isel(t=k)    # 1D array
Ek  = np.abs(np.fft.rfft(Ex))**2
```

---

## 8. Where this sits on the `PIC_NOTES.md` ladder

It is essentially **"Run 1.5"**: a driven electrostatic 1D run. Conceptually
two-stream-adjacent (1D, periodic, electrostatic, `B` stays boring) and it reuses
two example pgens almost verbatim. What makes it harder than Run 1 is *not* the
setup but three operational facts:

1. a real second species at `mᵢ/mₑ = 1836`;
2. the shot-noise budget forcing very large `ppc0` or a quiet start;
3. `10⁶`-step runs that put numerical heating and float precision front-and-centre.

None of that is GR (Run 4) or non-periodic-boundary machinery (Run 3).

### Run sequence — status

| # | run | params | status |
|---|---|---|---|
| 1 | undriven calibration | `A0 = 0`, `ppc0 ∈ {500,2000,8000}` | ✅ done — §0 |
| 2 | drive validation | `v_q/vₜₕ = 0.03`, short | ✅ done — `C = 0.1243` |
| 3–4 | amplitude ladder `run34` | 6 amplitudes + `Np` trio, `ωₚt = 3000` | ✅ done — 8 runs |
| — | **quiet start** implementation | `setup.loading = "quiet"`, 1V | ✅ done — §3b |
| R1 | **paper fast case** | `v_q/vₜₕ = 3e-2`, `N_pc = 200`, `ωₚt = 1e4` | ✅ **done** |
| R2 | **paper slow case** | `v_q/vₜₕ = 1e-3`, `N_pc = 2e5`, `ωₚt = 1e4` | ✅ **done** (48 h + 14 h restart) |
| H | LZ heating gate | undriven, `ppc0 = 2700`, `ωₚt = 2e4` | ✅ done — §0f |
| L1 | **LZ `v_q/vₜₕ = 1.0`** | `N_pc = 2.7e3`, 1 GPU | ▶ running |
| L2 | **LZ `v_q/vₜₕ = 0.1`** | `N_pc = 2.7e5`, 3 GPUs, ~81 h chained | ⬜ written, gated on L1 |
| 6 | 2D cross-check (optional) | one `A0` | ⬜ |
| 7 | long-term stability (optional) | `ωₚt → 8×10⁴` | ⬜ — needs multi-GPU |

Scripts: `run_paper_fast.sh`, `run_paper_slow.sh`, `run_paper_slow_resume.sh`,
`run_lz_heating.sh`, `run_lz1p0.sh`, `run_lz0p1.sh` in `pgens/dpdm/`.

### Results — the two resonance cases (paper Fig. 2)

Both at `ωₚt = 10⁴`, quiet start, 1V, single A100.

| | `N_pc` | `C` (theory 0.125) | peak `ε_E/ε_th` | at `ωₚt` | final `ΔKEₑ/ε_th` |
|---|---|---|---|---|---|
| fast `v_q/vₜₕ = 3e-2` | 200 | 0.1264 | 43.6 | 553 | **69.6** |
| slow `v_q/vₜₕ = 1e-3` | 2e5 | 0.1237 | 2.1 | 3694 | **2.4** |

`kTₑ(0)` comes out at `4.9963e-4` against the 1V target `5.0e-4`; the −0.075% is
**not** loading error but the relativistic definition `T00 = Σ w m γ`, whose next term
is `−(3/8)⟨u⁴⟩ = −0.075%`. The ions, 1836x heavier and far more non-relativistic,
land at `4.999998e-4` — which is the check that confirms it.

**The regime distinction lives in the final electron heating, not the field peak.**
Both cases overshoot `ε_thermal` in `ε_E` (43.6x vs 2.1x), but the lasting electron
heating differs by ~29x: `O(10²)` for fast, "not even changed by a factor of 2" for
slow — the paper's own phrasing for the LZ case. An earlier reading of the run34 data
claimed the slow regime "saturates *below* thermal"; at `N_pc = 2e5` with a quiet
start it does not — it peaks at 2.1x and then collapses.

**Saturation is delayed relative to run34.** With the noise floor 11 orders down at
`t = 0`, daughter waves grow from a far smaller seed, so saturation occurs later than
in the random-start runs — `ωₚt = 553` for the fast case rather than `~110`
(`t_sat + 1/ωₚᵢ`). Expected, and it means run34 timings are not the comparison.

## 9. Open risks — mostly closed

**Closed.**

- ~~Energy conservation over `4×10⁶` steps~~ — `ΔT00 = −½Δ⟨E²⟩`, `r = −0.996`; no grid
  heating. Confirmed again by §0f: `kTₑ` drifts `−0.016%` over 10⁶ steps with a quiet
  start.
- ~~Quiet start required and unimplemented~~ — built, validated, in production (§3b).
- ~~`ExternalFields` half-step phase~~ — `pusher_ctx.time` is `tⁿ`, correctly centred.
- ~~`ωₚ` calibration~~ — detuning scan, symmetric about `ω = 1`.
- ~~`Np` normalization in Eq. B4~~ — resolved from the paper source: `Nₚ` is the
  **total** macroparticle count. Eq. B4 is correct to 3.1x (§0c).
- ~~Filtering vs. daughter-wave fidelity~~ — `current_filters = 0` throughout; the
  daughters appear below `kλ_D = 1` as they should.

**Still open.**

- **Post-saturation cost.** §0g: per-step cost roughly doubles once the plasma goes
  turbulent, and the transition is not predictable from a short benchmark. Size
  walltimes with a 2x allowance and keep checkpoints on. This is now the main
  scheduling risk for every long run.
- **Collisions.** Entity (like SHARP) is collisionless; the paper treats the final
  ion-acoustic damping analytically via `τei`. Not a gap — it matches the paper's
  regime — but "persists forever" means "persists on collisionless timescales".
- **`ωₚt = 8×10⁴` long-term stability** (the paper's strongest claim) is untested. At
  `N_pc = 2e5` that is ~8x the slow case, so multi-GPU + chaining, and it is the one
  run where the post-saturation cost really bites.
- **No comparison against the published figures.** Every number here is internally
  consistent and matches the paper's *text*, but the figure PDFs in the source tarball
  have not been digitized. Peak `ε_E/ε_th` and its timing are the sharpest available
  discriminators.
- **2D cross-check** never run: 1+1D sufficiency is assumed on the paper's word.

## 10. File pointers

| What | Where |
|---|---|
| external-field hook example | `examples/external_fields/pgen.hpp`, `.toml` |
| external field applied to particle force | `src/kernels/pushers/sr.hpp:197-213` |
| pusher normalization `normalized_dt_half` | `src/kernels/pushers/sr.hpp:111` |
| Vay pusher | `src/kernels/pushers/sr.hpp:370` |
| streaming pgen (two-species Maxwellian pattern) | `pgens/streaming/pgen.hpp`, `twostream.toml` |
| density-per-pair splitting | `pgens/streaming/pgen.hpp:109-136` |
| uniform injector (random positions) | `src/kernels/injectors.hpp:83-289` (positions at `:202`) |
| co-located pair injection | `src/kernels/injectors.hpp:251-286` |
| `InjectUniformMaxwellians` | `src/archetypes/utils.h:47` |
| explicit-coordinate injection (`InjectGlobally`) | `examples/external_fields/pgen.hpp` `InitPrtls` |
| `T00` moment definition (`Σ w m γ`) | `src/kernels/particle_moments.hpp:214` |
| shape-function orders `S1..S11` | `src/kernels/particle_moments.hpp:118-138` |
| `StatsID` enum | `src/global/enums.h:286` |
| `FldsID` enum | `src/global/enums.h:248` |
| stats TOML syntax + `T` index note | `input.example.toml:509`, `:620-645` |
| `[scales]` unit derivation | `src/framework/parameters/parameters.cpp:47-78` |
| pgen concepts (`HasExternalFields`, ...) | `src/global/traits/pgen.h` |
| specialization registry (1D SRPIC Minkowski) | `src/framework/specialization_registry.h:41` |
| CMake `shape_order` / `deposit` handling | `CMakeLists.txt:38-95`, `cmake/defaults.cmake` |
| quiet-start kernel + `probit` | `pgens/dpdm/pgen.hpp` (`QuietInit_kernel`, `QuietNorm_kernel`) |
| quiet-start validation | `pgens/dpdm/validate_quiet.sh`, `check_quiet.py` |
| MPI decomposition-independence test | `pgens/dpdm/validate_mpi.sh` |
| production decks | `pgens/dpdm/run_paper_{fast,slow}.sh`, `run_lz{1p0,0p1}.sh` |
| per-kernel timing diagnostic | `pgens/dpdm/diag_slow.sh` (`blocking_timers`) |
| run34 extraction + plots | `pgens/dpdm/run34/extract_run34.py`, `run34_plots.ipynb`, `PLOTS.md` |
| paper source (all parameters) | `arXiv-2510.13956v1.tar.gz` → `main.tex` |
