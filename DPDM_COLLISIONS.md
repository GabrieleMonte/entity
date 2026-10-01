# Collisional effects in DPDM resonant conversion

The case for adding a Coulomb collision operator to the DPDM problem generator,
and a plan for building one in Entity.

Companion to `DPDM_PIC_NOTES.md` (the reproduction plan for Hook, Huang & Shalaby,
arXiv:2510.13956 — "HHS" throughout). That file establishes that we reproduce the
paper. This file is about the one piece of physics the paper leaves out.

> **Status (2026-09-30).** Analysis only. No collision code written yet. Every
> number in §2 and §3 is *measured* or *computed*, and reproducible with the two
> scripts in `pgens/dpdm/collisions/` — see §8. §4 is argued, not measured, and
> flags its own uncertainties. §6 is a plan, not a result.
>
> **The short version.** HHS shut off resonant conversion with a saturated state
> whose stability rests on one unexamined clause: *"the damping rate on electrons
> would be polynomial suppressed **if the electrons were thermalized with
> themselves**"* (`main.tex:203`). That clause is false in their simulations, and
> it cannot be otherwise: the runs are 1+1D, and **electron–electron collisions
> are mathematically impossible in 1V**. Their longest run covers **0.8% of one
> e–e collision time**; the physical resonance lasts **~2×10⁴ of them**.

---

## 1. What the paper actually argues

Worth getting exactly right, because the opening is narrower and sharper than
"they ignored collisions".

**They do not dismiss electron–ion collisions.** What `main.tex:203` calls
exponentially suppressed is *ion Landau damping of the ion-acoustic waves*,
`exp[-(c_s/v_th^i)²/2] ≈ exp[-T_e/2T_i]` — a wave-damping statement, not a
collision statement. On actual e–i collisions they do the opposite of dismissing:
they write down `ν_ei` (`main.tex:207`), compute the energy-exchange time
`τ_ei ≈ 10⁴ s` at recombination (`main.tex:215`), note that e–i collisions reduce
`T_e/T_i` and sustain the steady state (`main.tex:628`), and close with:

> *"Dedicated simulations of a collisional plasma might be useful to understand
> this final state."* (`main.tex:216`)

That is an open invitation, and it is the work proposed here.

**The actual gap is electron–electron collisions.** The shutoff argument depends
on the electrons being self-thermalized, and:

- `ν_ee` appears **zero times** in the paper. `ν_ei` and `ν_eγ` are both tabulated
  and used throughout the non-resonant-heating appendix; the electron
  self-collision rate — the one that decides whether the clause holds — is absent.
- The phrase *"thermalized with themselves"* occurs **once**, as a conditional,
  and is never revisited.
- The simulations are **1+1D** (`main.tex:133`, `main.tex:485`): only the `(x,v_x)`
  subspace. In 1V an elastic collision between two identical particles has only
  two solutions — identity, or velocity exchange, which for indistinguishable
  particles is the same state. **e–e collisions are not weak in a 1V code; they
  do not exist.** The single process their argument relies on is the one their
  setup structurally cannot contain.

## 2. The conditional is false — measured

From `fast_s` (decohering loader, `A0 = 9.48683e-4`, `Nx = 1000`, `ppc0 = 200`),
last 6 phase-space frames, `ω_p t ≈ 9150–10⁴`, integrated over `x`, compared
against a Maxwellian of the **same variance**. Reproduce with `fdist_diag.py`.

### Electrons: excess kurtosis **+6.2**, core deviates **57%**

| v/v_th | counts | stat err | f / Maxwellian |
|---:|---:|---:|---:|
| 0.26 | 2033 | 2.2% | **1.68** |
| 0.52 | 1055 | 3.1% | 0.97 |
| 0.99 | 301 | 5.8% | **0.39** |
| 1.51 | 136 | 8.6% | **0.34** |
| 1.98 | 80 | 11.2% | 0.46 |
| 2.51 | 44 | 15.1% | 0.82 |
| 2.98 | 38 | 16.2% | **2.57** |
| 3.50 | 27 | 19.2% | **9.93** |
| 4.02 | 19 | 22.9% | **50.0** |
| 4.50 | 10 | 31.6% | **195.7** |

Sharply peaked core, depleted shoulder at 1–2 `v_th`, suprathermal tail enhanced
by 1–2 orders of magnitude. Per-bin errors are large in the tail, but the trend
spans a factor of 200 against 20–30% noise — it is real. A core-plus-halo
structure, frozen in by the absence of any entropy-producing process.

**Do not read a bump-on-tail into this.** Smoothing gives 8 positive-slope bins
out of 80 in `0.3 < w < 4.5`, which is consistent with the per-bin noise at those
counts. The *systematic* shape is robust; individual bins are not.

### Ions: excess kurtosis **+1.6**, and we cannot yet see the velocity that matters

`f_i/Maxwellian` rises 0.98 → 0.79 → 0.64 → 0.88 → 1.80 → 3.68 → **17.6** at
1, 1.5, 2, 2.5, 3, 3.5, 4 `v_th^i`.

**Known gap.** `c_s/v_th^i = 7.31` but the extractor's ion histogram clips at
5.3 `v_th^i`, so `f_i(c_s)` — the quantity that actually sets ion Landau damping
of the IAW — **cannot currently be measured**. The `exp(-26)` suppression HHS
quote is a *Maxwellian* result, and the ion distribution is visibly not
Maxwellian in the direction that matters. Widening the histogram range is a
prerequisite (§6, Gate 0).

### `T_e/T_i` never saturates

| ω_p t | 100 | 300 | 1000 | 3000 | 5000 | 9000 |
|---:|---:|---:|---:|---:|---:|---:|
| T_e/T_i | 2.5 | 12.7 | 25.4 | 37.8 | 43.3 | **52.8** |

HHS quote "electrons about 30 times hotter than the ions" as the steady state.
Ours passes through 30 near `ω_p t ≈ 1500` and is **still rising monotonically**
at the end of the run. Since ion Landau damping goes as `exp(-T_e/2T_i)`, this is
a **self-reinforcing lock**: the hotter the electrons get relative to the ions,
the deader the only dissipation channel becomes, the more the energy stays stuck.
Nothing in a collisionless 1V plasma ever breaks the loop. "30" is not a steady
state — it is where the run stopped.

### Absorption never actually stops

Residual `d(ΔKE_e)/dt` in `fast_s`, code units per `ω_p t`:

| window | rate |
|---|---|
| [1000, 5000] | 3.62e-6 |
| [5000, 9500] | 1.85e-6 |

Decaying slowly, nonzero, and consistent with HHS's own corrected curve
(1.91e-6 — see the factor-of-2 finding in `DPDM_PIC_NOTES.md`). The system is
not frozen; it is leaking.

## 3. The decisive argument — they simulate 0.8% of one collision time

Rates from standard cosmology, **calibrated against the paper's own number**:
`cosmo_rates.py` gives `τ_ei(z=1100) = 1.24e4 s` against their quoted
"about 10⁴ s" (`main.tex:215`). Normalization confirmed.

The key structural result: **`ν_ee/ω_p` is epoch-independent in the early
universe.** Since `ν ∝ n^{1/2} T^{-3/2}` with `n ∝ a⁻³` and `T ∝ a⁻¹`, the ratio
is flat from `z = 10⁶` down to recombination:

```
nu_ee / omega_p  ~=  1.0e-7        (equivalently  n_e * lambda_D^3 ~ 3e6)
```

The full hierarchy, in `ω_p t`:

| | ω_p t | in e–e collision times |
|---|---:|---:|
| ion response `1/ω_p^i` | 43 | — |
| burst saturation | ~300–500 | — |
| our runs | 1e4 | 1e-3 |
| **HHS longest run** | **8e4** | **0.008** |
| **e–e collision time** | **1e7** | 1 |
| **e–i energy exchange `τ_ei`** | **9e9** | 650 |
| **resonance duration `ε/H`** (ε=1e-8) | **2e11** | **~2e4** |

**HHS's longest simulation covers 0.8% of a single electron–electron collision
time. The physical resonance lasts roughly twenty thousand of them.** Both
collision timescales fall strictly *inside* the resonance window and strictly
*outside* the simulated window.

The conclusion that resonant conversion deposits only O(1–100) × `ε_th` — which
is what removes the DPDM constraints in their Fig. 5 — is an extrapolation of a
collisionless frozen state across six orders of magnitude in time, through a
regime where the thing holding it frozen is guaranteed to have been erased ~10⁴
times over.

### Two corollaries

- **`ν_ii ≈ 9 ν_ee` in the saturated state.** `ν_ii/ν_ee = (m_e/m_i)^{1/2}
  (T_e/T_i)^{3/2} = 9.0` at `T_e/T_i = 52.8`. **Ions are the most collisional
  species** precisely in the state the shutoff requires. Any collisional
  treatment must include i–i, not just e–e.
- **Their collision-relevance remark does not survive the arithmetic.**
  `main.tex:203` says `ω_p t = 8e4` "is already approaching when collisions
  between particles shall become relevant". By this calculation it is short by
  ~125× for e–e and ~10⁵× for e–i energy exchange. ⚠ **Verify independently
  before this goes into any correspondence with the authors** — but the
  normalization is pinned by their own `τ_ei`.
- **Compton is faster still at the earliest times.** `ν_eγ/ω_p ∝ a^{-5/2}`,
  reaching ~1e-6 at `z = 10⁶` — ten times `ν_ee/ω_p`. Not modelled here, but it
  only strengthens the case at high `z`.

## 4. Which way does the effect go? (argued, not measured)

Ranked by defensibility. **The sign is not guaranteed.**

**M1 — Longitudinal energy drain by isotropization. Strongest; a direct
consequence of their 1V choice.** In 1D3V, collisions move energy from the heated
parallel degree of freedom into two cold transverse ones at rate `~ν`. The
parallel temperature sets the relativistic detuning, the nonlinearity parameter
`v_q/v_th∥`, and `T_e/T_i`. Full isotropization drops `v_th∥` by `√3`
(0.29c → 0.17c) *at fixed total energy*. All three push back toward absorbing.
⚠ **Caveat on the record:** in real 3D the turbulence isotropizes on its own, so
part of M1 is a 1D artifact. But the constraint comes from a 1D1V run, so testing
its sensitivity to exactly this is legitimate.

**M2 — e–i energy exchange unlocks ion Landau damping. Strong; their own admitted
mechanism.** Driving `T_e/T_i` from 53 toward 1 changes `exp(-T_e/2T_i)` by ~10¹¹,
draining the ion-acoustic energy that constitutes the blocking turbulence.
`ω_p τ_ei = 9e9 < 2e11`, so it completes ~30× during the resonance.

**M3 — Restoring `∂f/∂v` where the Langmuir daughters live. Ambiguous sign —
measure, do not predict.** We see 50–200× Maxwellian at 3–4 `v_th`, where the
`kλ_D ~ 0.3` daughters resonate. Collisions deplete that tail, which *reduces*
Langmuir damping (more turbulence survives → hurts), but also recycles the energy
into the bulk and resets the shape (→ helps).

**M4 — Relativistic detuning. Collisions cannot undo it.** Collisions conserve
energy, so they cannot cool the plasma back onto resonance. They change only the
*shape* at fixed energy, and `ω_p,eff² = ω_p²⟨γ⁻³⟩` is shape-dependent at the
few-% level. Small, sign unclear, and the physical resonance is `~ε = 1e-8` wide,
which makes any detuning argument delicate. **Do not build the case on this.**

**Verdict.** The mechanism-level case is suggestive but not airtight. The
timescale case (§3) is airtight and sufficient on its own. **A null result is
still worth having** — it would convert HHS's assumption into a validated one.

## 5. Why 1V is a showstopper

`setup.one_v = true` (`pgens/dpdm/pgen.hpp:316`, applied at `:247`) zeroes `ux2`
and `ux3` to match SHARP's 1V Maxwellian. Under it:

- **e–e collisions:** impossible (§1).
- **i–i collisions:** impossible, same argument.
- **e–i collisions:** exist (unequal masses give a nontrivial 1D elastic
  solution) but are *deterministic* — a single outcome, not a Coulomb scattering
  distribution. Not physical either.

**1D3V is a hard prerequisite.** The flag flip is trivial; the consequences are
not (§6, Gate 0).

## 6. Implementation plan

### What Entity already has

`grep -ril collision src/` returns **nothing** — no collision module exists. But
the two hard prerequisites do:

| need | where |
|---|---|
| per-cell particle binning | `Particles::SortSpatially`, `src/framework/containers/particles_sort.cpp:197` (`Kokkos::BinSort` over cell index) |
| called every step in srpic | `src/engines/srpic/srpic.hpp:185` |
| a hook right after the sort | `CustomPostStep(step, time, domain)`, declared `src/global/traits/pgen.h:109`, called `src/engines/engine.hpp:276` |
| thread-safe RNG | `domain.random_pool()` |

Ordering is `step_forward` (ends with `SortParticles`) → `CustomPostStep`. So the
whole operator can live in `pgens/dpdm/pgen.hpp` with **zero changes to Entity
core**, which keeps the `dpdm` branch rebaseable. Upstreaming to `src/kernels/`
is a later decision.

### Gate 0 — 1D3V baseline (blocking)

1. `setup.one_v = false`; re-run `fast` **collisionless** as the control.
   Transverse momenta are passive in a 1D electrostatic run *except* through `γ`,
   and at post-burst `v_th ~ 0.3c` that coupling is real — so the 3V control is
   genuinely a different run and must be measured, not assumed.
2. Widen the phase-space histogram `u`-range in `extract_paper.py` so ions reach
   `≥10 v_th^i` (currently clips at 5.3; `c_s` is at 7.3 — see §2).
3. Check `n_dof` bookkeeping end-to-end: `ε_th` is 3× larger in 3V, and
   `read_stats(tag, n_dof, ...)` already takes it as a parameter.

Cost ~3 h, 1 GPU. ⚠ **This moves the baseline off the paper's exact 1V setup.**
Every collisional result will be against a 3V control, so the 3V collisionless
run must be documented as a result in its own right, not as scaffolding.

### Gate 1 — numerical collisionality floor (blocking)

A physical operator is meaningless if PIC noise already relaxes `f` faster. Load
a deliberately non-Maxwellian `f_e` (two-temperature, or the measured core+halo),
run **undriven** with collisions **off**, fit the relaxation rate → `ν_num`.
**Require `ν_phys ≳ 10 ν_num` for every production run.**

Encouraging prior: the `npc` knob run (10× particles, same `dx`) left the late
absorption rate unchanged (4.40e-6 vs 3.82e-6 base, inside realization scatter),
which says the current late-time behaviour is *not* numerical-collision driven.
Cost <1 h.

### Step 2 — the operator

**Takizuka–Abe (1977) pairwise Monte Carlo, with the Perez et al. (2012,
Phys. Plasmas 19:083104) relativistic extension.** Entity is `srpic`/Vay and
post-burst `v_th ~ 0.3c`, so the relativistic form is required; it reduces to
plain TA in the cold limit.

**The decisive design choice: TA conserves momentum and energy exactly, pair by
pair.** The entire measurement is an energy budget at the 1e-6 level (§2) — a
Langevin/Fokker–Planck operator's energy non-conservation would swamp the signal
outright. Do not compromise on this.

- Hook: `CustomPostStep`, with `spatial_sorting_interval = 1`.
- Per-cell ranges: after the sort `i1` is monotonic, so recompute cell offsets
  with a trivial 1D scan in the pgen rather than exposing `BinSort`'s bin offsets
  from core.
- Pairs: **all three** — e–e, i–i (recall `ν_ii ≈ 9 ν_ee`, §3), e–i. Random
  shuffle within cell, TA pairing, TA's 3-particle trick for odd counts.
- Assert uniform weights, or use Perez's unequal-weight pairing.
- Knob: `setup.nu_ee_norm` = `ν_ee/ω_p` at the **initial** temperature,
  implemented as a multiplier on `lnΛ` so the physical `T^{-3/2}` scaling is
  preserved as the plasma heats.
- Stability: `ν dt = 1e-4` at `ν/ω_p = 1e-2`, `dt = 0.01`. Comfortable.

### Step 3 — the physics runs

`Nx = 2000` (our production rule for post-burst energies), to `ω_p t = 1e4`,
scanning `ν_ee/ω_p ∈ {0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2}`. The top of that range
gives ~100 collision times in one run — a genuinely collisional steady state.

**The key test is not "does it change".** It is whether the post-burst deposition
**collapses onto a function of `ν_ee·t` alone**. If it does, extrapolation to the
physical `ν t ~ 2e4` rests on a measured scaling law instead of an assumption —
which is the whole point, and exactly what the paper lacks.

Observables: post-burst `d(ΔKE_tot)/dt` vs `ν`; `T_e/T_i(t)` — does the monotonic
climb break?; excess kurtosis and tail ratio vs `t`; `f_i(c_s)`.

**Run one variant with collisions switched on only at `ω_p t = 1000`,** after the
burst. At `ν/ω_p = 1e-2` collisions are not a perturbation *during* the burst, so
this cleanly isolates the effect on the saturated state without contaminating the
physics that produced it.

Cost, from the knob-run accounting (`k_nx` at `Nx=2000` was 1h26, `k_npc` 1h53):
~4–5 h per run with collisions, **~30 GPU-hours for the full scan**. `normal`
queue.

## 7. Risks

- **1D is 1D.** 1D Langmuir/ion-acoustic turbulence is not 3D turbulence, and M1
  is partly a 1D artifact (§4).
- **Artificially enhanced collisionality is a model, not the real plasma.** The
  scan-and-collapse (§6 Step 3) is what makes it an argument rather than a
  demonstration. We cannot reach `ν/ω_p = 1e-7` directly — that is the whole
  reason the paper's gap exists.
- **The sign is genuinely not guaranteed** (§4).
- **Gate 0 moves the baseline** off the paper's 1V setup (§6).

## 8. Reproducing the numbers

```sh
module load python/3.12.11
python3 pgens/dpdm/collisions/cosmo_rates.py                       # Section 3
python3 pgens/dpdm/collisions/fdist_diag.py \
        $SCRATCH/dpdm_npy/fast_s_phase.npz 6                       # Section 2
```

`fast_s_phase.npz` is produced by `pgens/dpdm/paper/extract_paper.sh`.
`cosmo_rates.py` self-checks against the paper's quoted `τ_ei`.

## 9. Open questions

1. `f_i(c_s)` is unmeasured (histogram clipping, §2). Until it is, the claim that
   ion Landau damping is *not* as dead as `exp(-26)` suggests is extrapolation.
2. Does the 1V→3V change alone alter the saturated state? Unknown until Gate 0.
3. Is the residual absorption (§2) collisionally limited, or something else? Our
   knob scan found it `dx`-dependent with **no identified mechanism**
   (see `DPDM_PIC_NOTES.md`); that remains open and could confound the collision
   scan. Gate 1 partially addresses it.
4. Does `ν_eγ` (Compton) need to be in the model for the early-`z` case?
5. Should this be upstreamed to `src/kernels/` as a general Entity feature?
