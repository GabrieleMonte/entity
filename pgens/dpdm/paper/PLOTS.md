# paper runs — plots to make, and what each one claims

These six runs are the paper-faithful production set. They are **not** a subset of
`run34/`: that survey is random-start at a uniform `ppc0 = 2000` out to `ωₚt = 3000`.
These are quiet-start, at the paper's own `N_pc`, out to `ωₚt = 10⁴`–`2×10⁴`, and
include the swept-frequency family that no resonance run can produce.

| tag | family | `v_q/vₜₕᵉ` | `A0` | `N_pc` | `ωₚt` | loading |
|---|---|---|---|---|---|---|
| `fast` | resonance | 3e-2 | 9.48683e-4 | 2e2 | 1e4 | quiet |
| `slow` | resonance | 1e-3 | 3.16228e-5 | 2e5 | 1e4 | quiet |
| `lz1p0` | landau_zener | 1.0 | 3.16228e-2 | 2.7e3 | 2e4 | quiet |
| `heat` | undriven | — | 0 | 2.7e3 | 2e4 | quiet |
| `undriven_q` | undriven | — | 0 | 64 | 400 | **quiet** |
| `undriven_r` | undriven | — | 0 | 64 | 400 | **random** |

All: `Nx=1000`, `L=40 c/ωₚ`, `dt=0.02`, `mᵢ/mₑ=1836`, `kTₑ=kTᵢ=10⁻³ mₑc²`, 1V.

Constants: `vₜₕᵉ = 0.0316`, `ε_th(1V) = 5e-4`, `1/ωₚᵢ = 42.8`, `kλ_D = 1` at `k = 31.6`,
`ε_noise = 2.06e-4/N_pc`. Every per-run constant is in `meta.json`.

---

## Plot A — energy budget, paper Fig. 2. **The central figure.**

`fast_stats.npz`, `slow_stats.npz`. Two panels, log–log vs `ωₚt`.

Curves per panel: `eps_E`, `dKE_e`, `dKE_i`, `eps_tot`, the line `(A0²/8)(ωₚt)²`, and
horizontals at `eps_thermal_1V = 5e-4` and `eps_noise`.

**Claim:** `t²` growth, then saturation, in both regimes — but the *outcomes* differ.
Measured: `fast` peaks at `eps_E/ε_th = 43.6` at `ωₚt = 553` and ends with
`dKE_e/ε_th = 69.6`; `slow` peaks at `2.1` at `ωₚt = 3694` and ends at `2.4`.

> **The discriminator is the final electron heating, not the peak field energy.**
> Both overshoot `ε_th` in the *field*; the lasting heating differs by ~29×. That is
> the paper's own phrasing for the slow case ("the kinetic energy of the electrons is
> not even changed by a factor of 2").

> Draw the linear line against **`eps_tot`**, never `eps_E` — factor 2. Thermal
> reference is the 1V `5e-4`, not `T00_1 − 1 = 1.5e-3` which is 3V.

This supersedes `run34` Plot 1: same shape, but at the paper's parameters and with a
quiet start. Note the **saturation is later** than in the random-start run34 scans —
`ωₚt = 553` rather than ~110 — because daughter waves grow from a floor eleven orders
lower. Random-start timings are not a valid comparison.

## Plot B — the Landau–Zener transition (paper Fig. 3). **New physics.**

`lz1p0_stats.npz`, `lz1p0_temps.npz`, `meta.json` → `omega0`, `domega_dt`, `t_cross`.

`eps_E`, `dKE_e`, `eps_tot` vs `ωₚt`, with `ω(t)/ωₚ = 0.8 + 2e-5·t` on a twin axis.
Mark two times: `ωₚt = 5000` (`ω = 0.9ωₚ`) and `ωₚt = 10⁴` (`ω = ωₚ`, **the crossing**).
Overlay `e_ext` from the stats file to confirm the swept phase is what the pusher
actually used.

**Claim:** conversion is killed by nonlinearity *before* the crossing. The paper's
"nonlinearity kills the resonant transfer around `ωₚt = 5000`" refers to `ω = 0.9ωₚ`
— off-resonant excitation, not the crossing itself. `dKE_e` should not change by a
factor of 2 across the whole sweep.

> Unlike the resonance family, the pre-crossing drive is a **bounded** forced
> oscillation of amplitude `A0/(1−ω²/ωₚ²) = 2.78·A0` at `ω = 0.8ωₚ`. There is no time
> at which a growing signal outruns the noise, which is why `N_pc` is set by
> `t̄_noise ≤ 0.02` rather than by a saturation argument.

## Plot C — quiet vs random loading. **The method gate.**

`undriven_q_*` and `undriven_r_*` overlaid. Identical decks; `loading` is the *only*
difference. Two panels:

1. `eps_E(t)` from `_stats.npz`, log-y, both runs, with `eps_noise = 3.22e-6` marked.
2. `dn_over_n` from `_fields.npz`: rms over `x` vs `t`, both runs.

```python
rms = f["dn_over_n"].std(axis=-1)      # (nt,)
```

**Claim:** the quiet lattice starts at round-off and *stays* there under free
streaming, while the random start sits at the shot-noise floor from `t = 0`. Verified
in this pipeline: rms `δn/n` is **exactly 0.0** for the quiet run at both `t = 0` and
`t = 400`, against `0.129 ≈ 1/√64 = 0.125` for random — i.e. the random run is at
textbook shot noise and the quiet run has none. `eps_E` differs by ~13 orders.

> **Do not use `Charge` or `Ex(t=0)` as the discriminator.** Both are zero for *both*
> loaders — Entity initialises `E` to zero identically and both co-locate the species.
> The channel that matters is the **total** density `n_e + n_i`, which is what
> `dn_over_n` holds, because `ωₚ ∝ √n`.

## Plot D — what the quiet start does *not* buy. **The sizing result.**

`heat_stats.npz` (undriven, `N_pc = 2.7e3`, `ωₚt = 2×10⁴`, 10⁶ steps). Two panels:

1. `eps_E(t)`, log-y: starts at `6.7e-23`, climbs, plateaus at `4.8e-8` by
   `ωₚt ≈ 5000`. Mark the random-start floor at the same `N_pc`, `7.6e-8`.
2. `Te_rel`, `Ti_rel` — linear y, zoomed. Drift over the full run is
   `−0.0162%` / `+0.0084%`.

**Claim, panel 1:** the ordered state **relaxes to ~62% of the thermal shot-noise
floor**. A quiet start buys a clean `t = 0` and nothing more — it does **not** permit
a reduction in `N_pc`, which must still be sized from Eq. B4. This is the paper's own
"it takes `tωₚ ~ 3–4` for the driven energy to be converted into de-coherent high-`k`
noise", seen from the undriven side.

**Claim, panel 2:** an ordered lattice does not drift. This is the gate that makes the
Landau–Zener runs meaningful — a slow `Tₑ` creep would shift `ωₚ` and masquerade as
the sweep crossing resonance.

## Plot E — `E(k,t)` at paper parameters. **The mechanism.**

`fast_fields.npz`, `slow_fields.npz`: `PE`, `Pn1`, `Pn2`, `k`, `w_k`.

`(k, ωₚt)` heatmap of `log10 PE`, plus the `k=0` / `k≠0` split on Plot A's axis:

```python
P_k0  = w_k[0] * PE[:, 0]
P_knz = (w_k[1:] * PE[:, 1:]).sum(axis=1)     # sum(w_k*PE) == mean(Ex^2), exactly
```

**Claim:** `k≠0` power starts rising exactly where `eps_tot` leaves the linear line in
Plot A — modulational daughters in `PE`, the ion-acoustic branch at low `k` in `Pn2`.
Check the daughters live below `kλ_D ~ 1` (`k < 31.6`); action at the grid scale would
be numerical, not physical.

Same construction as `run34` Plot 2, but at the paper's `N_pc` and with the quiet
start, so the *onset* is trustworthy rather than seeded by initial noise.

## Plot F — phase space. **The one to animate.**

`<tag>_phase.npz`: `H_e`, `H_i` of shape `(60, 240, 240)`, with `x_edges`,
`u_edges_e`, `u_edges_i`, `t`. 60 frames spanning the run, fixed axes.

```python
plt.pcolormesh(ph["x_edges"], ph["u_edges_e"], ph["H_e"][i].T, ...)
```

For `fast`: a flat thermal band, then the driven Langmuir wave folding into vortices
at saturation, then a hot disordered distribution. For `lz1p0`: the same, arriving as
`ω(t)` approaches `ωₚ`. Pair frames with the `Ex(x)` and `N_1+N_2` snapshots in
`_fields.npz` at matching times to show the density cavitation directly.

**Claim:** density cavitation makes `n_tot(x)` — hence `ωₚ(x)` — non-uniform, which is
what detunes the resonance and shuts conversion off. This is the visual version of the
argument Plots A and E make quantitatively.

---

## What counts as a successful reproduction

1. `eps_tot` tracks `(A0²/8)(ωₚt)²` before saturation in both cases (Plot A).
2. `fast` overshoots `ε_th` and ends `O(10²)` hotter; `slow` barely heats (Plot A).
3. Landau–Zener transfer dies before the crossing at `ωₚt = 10⁴` (Plot B).
4. Quiet start is at round-off and random is at `1/√N_pc` (Plot C).
5. The ordered state relaxes to the shot-noise floor by `ωₚt ≈ 5000` (Plot D).
6. `k≠0` power rises at the departure point, below `kλ_D ~ 1` (Plot E).

1–3 are the paper's core argument. 4–5 validate the method. 6 is the mechanism.

## Known caveats

- **`lz0p1` is missing.** The low-amplitude Landau–Zener case (`v_q/vₜₕᵉ = 0.1`,
  `N_pc = 2.7e5`, 3 A100s, ~81 h chained) is written but gated. Plot B is one
  amplitude, not the family.
- **No comparison against the published figures.** Everything here is checked against
  the paper's *text*; the figure PDFs in the arXiv source have not been digitised.
  Peak `eps_E/ε_th` and its timing are the sharpest available discriminators.
- **`ωₚt = 10⁴`, not `8×10⁴`.** The long-term-stability claim is untested.
- **Collisionless**, like SHARP. "Persists indefinitely" means "on collisionless
  timescales"; the paper treats ion-acoustic damping via `τ_ei` analytically.
- **1V** (`one_v = true`, `uy = uz = 0`). Matches SHARP exactly, but it means every
  thermal reference is the 1V `5e-4` — see the trap notes in `extract_paper.py`.
