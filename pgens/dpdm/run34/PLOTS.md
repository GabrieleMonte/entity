# run34 — plots to make, and what each one claims

run34 is the **resonance family only** (`ω = ωₚ` fixed). Paper Fig. 3 (Landau–Zener)
is *not* reproducible from these runs — that needs the swept-frequency family,
step 5 of the plan.

Eight runs, all `Nx=1000`, `L=40 c/ωₚ`, `dt=0.02`, `mᵢ/mₑ=1836`, `kTₑ=kTᵢ=10⁻³ mₑc²`:

| tag | `v_q/vₜₕᵉ` | `A0` | `ppc0` | `ωₚt` | `t_sat` |
|---|---|---|---|---|---|
| `r1e-3` | 1e-3 | 3.162e-5 | 2000 | 6000 | 1999 |
| `r3e-3` | 3e-3 | 9.487e-5 | 2000 | 3000 | 666 |
| `r1e-2` | 1e-2 | 3.162e-4 | 2000 | 3000 | 200 |
| `r3e-2` | 3e-2 | 9.487e-4 | 2000 | 3000 | 67 |
| `r1e-1` | 1e-1 | 3.162e-3 | 2000 | 3000 | 20 |
| `r3e-1` | 3e-1 | 9.487e-3 | 2000 | 3000 | 6.7 |
| `conv500` | 3e-3 | 9.487e-5 | 500 | 3000 | 666 |
| `conv8000` | 3e-3 | 9.487e-5 | 8000 | 3000 | 666 |

Constants: `vₜₕᵉ = 0.0316`, `ε_thermal(1V) = 5e-4`, `1/ωₚᵢ = 42.8`, `kλ_D = 1` at `k = 31.6`.

---

## Plot 1 — energy budget (paper Fig. 2). **The central figure.**

`<tag>_stats.npz`. Two panels: slow `r3e-3`, fast `r1e-1`. Log-log vs `ωₚt`.

Curves: `eps_E`, `dKE_e`, `dKE_i`, `eps_tot`, the line `(A0²/8)(ωₚt)²`, and a
horizontal `eps_thermal_1V = 5e-4`.

**Claim:** `t²` growth, then saturation. Slow saturates *below* the thermal line;
fast *overshoots* it and halts roughly `1/ωₚᵢ = 43` after `t_sat`.

> Draw the linear line against **`eps_tot`**, never `eps_E` — factor 2. And the
> thermal reference is the 1V `5e-4`, not `T00_1 − 1 = 1.5e-3` which is 3V.

## Plot 2 — `E(k,t)` and `n(k,t)` (paper Fig. 8, lower). **The mechanism.**

`<tag>_fields.npz`: `PE`, `Pn1`, `Pn2`, `k`, `w_k`.

`(k, ωₚt)` heatmap of `log10 PE`, plus the `k=0` vs `k≠0` split on Plot 1's axis:

```python
P_k0  = w_k[0] * PE[:, 0]
P_knz = (w_k[1:] * PE[:, 1:]).sum(axis=1)     # sum(w_k*PE) == mean(Ex^2)
```

**Claim:** `k≠0` power starts rising exactly where `eps_tot` leaves the linear
line in Plot 1 — modulational daughters in `PE`, the ion-acoustic branch at low
`k` in `Pn2`. Without this you've shown *that* it saturates, not *why*. Check the
daughters live below `kλ_D ~ 1` (`k < 31.6`); if the action is at the grid scale
it's numerical, not physical.

## Plot 3 — `Np` convergence. **The credibility gate.**

`conv500` / `r3e-3` / `conv8000` (`ppc0` = 500 / 2000 / 8000, same `A0`) overlaid
on Plot 1. Mark `eps_noise = 2.06e-4/ppc0` per run.

**Claim:** saturation is physical, not shot noise. If the three overlay, the quiet
start is unnecessary — which closes **step 4** of the plan and confirms §3b of the
notes against the paper's Eq. B4 estimate.

## Plot 4 — regime transition. **What the six-point ladder bought.**

`summary.npz`. `eps_E_max_norm` vs `vq_over_vthe`, log-log, with `boundary`
(`(mₑ/mₚ)^½/2 = 0.0117`) marked.

**Claim:** the slow/fast dichotomy is a genuine threshold, not two anecdotes —
points below the boundary stay under 1, points above overshoot.

## Plot 5 — `Tₑ/Tᵢ(t)` (paper Fig. 7).

`<tag>_temps.npz`: `Te`, `Ti`. Fast-growth runs.

**Claim:** `Tₑ/Tᵢ → O(30)`, electrons `O(10²)` times hotter than initial. Both are
**parallel** temperatures from the `ux` variance (§6.5) — a 3V temperature would
be wrong by 3.

## Plot 6 — real-space structure (paper Fig. 8, upper).

`<tag>_fields.npz`: `N_1`, `N_2`, `Rho`, `Ex` vs `x` at a few times bracketing
saturation.

**Claim:** density cavitation makes `n_tot(x)` — and hence `ωₚ(x)` — non-uniform,
which is what detunes the resonance and shuts conversion off.

---

## What would count as a successful reproduction

1. `eps_tot` tracks `(A0²/8)(ωₚt)²` before saturation, all amplitudes (Plot 1).
2. Slow saturates below `5e-4`; fast overshoots then halts (Plots 1, 4).
3. `k≠0` power rises at the departure point (Plot 2).
4. The three `ppc0` overlay (Plot 3).
5. `Tₑ/Tᵢ` reaches tens in the fast runs (Plot 5).

1–4 are the paper's core argument. 5 is corroboration.

## Known caveats

- **Collisionless.** Like SHARP. The paper's "persists forever" is really
  "persists on collisionless timescales" — ion-acoustic damping via `τei` is
  treated analytically in the paper, not simulated.
- **`ωₚt = 3000`, not `8×10⁴`.** The long-term-stability claim (run 7) is not
  tested here.
- **3V vs 1V** everywhere temperatures or thermal energies appear.
