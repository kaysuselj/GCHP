# Scope: an exact (discrete) adjoint of GCHP 14.8.0 tracer advection

Written 2026-10-16. Status: not started. Decision 2026-10-16: start with option B (linear-scheme transpose), development on the laptop from the week of 2026-10-19.

## 1. Why

The current advection adjoint (`ADV_ADJ_METHOD: reverse_wind_pure`, `ADV_ADJ_WEIGHT: airmass`) runs the
**forward** FV3 scheme with reversed winds on λ / (air mass). It is a continuous-adjoint approximation,
not the transpose of the forward operator. What the tests show:

| test | current result |
|---|---|
| P2a/P2d (8 m/s, face edge/corner) and P2c (tilted rotation) | errors concentrated at edges and corners |
| T2 blocks, real winds, 5 days | 4.4 / 5.5 % |
| P2e blocks, real winds, 1 month (airmass) | 4–13 % |
| OSSE FD check, 12 cells | slope 0.994, weighted error 1.07 % |
| linear-scheme dot-product test | error/‖λ‖ 27–190 % within 2 cells of face edges/corners |

The causes are known (agent bug hunt, 2026-10-15):
- The PPM edge formulas (`al(0..2)` in `xppm`/`yppm`) are not reflection-symmetric.
- The limiter and the mass fixer are nonlinear.
- `copy_corners` has a small transpose defect at the faces 0/2/4 and 1/3/5 corners.

None of these can be fixed inside the reversed-wind approach. Only a discrete adjoint removes them. With a
discrete adjoint, the FD agreement should reach ~1e-3 or better. The remaining limit is limiter switches
flipping under the FD perturbation. Agreement would also hold at C48/C90, where the edge cells are a
larger share of the observed area.

The forward model is not touched. The new adjoint would be a third option next to `reverse_wind_pure`.

## 2. What the forward operator is (14.8.0, our configuration)

One advection step, `offline_tracer_advection` (`fvdycore/model/fv_tracer2d.F90:846`), CO2 only (nq = 1):

1. **Halo setup of CX/CY/MFX/MFY** (`mpp_get_boundary`, `mpp_update_domains`). These are winds, not
   active variables, so no adjoint is needed: the fluxes are recomputed from met in the adjoint, as now.
2. **Horizontal transport** `tracer_2d` (`z_tracer` false):
   - sub-cycling `nsplt = int(1 + max Courant)`, a global max (discrete, not differentiated);
   - per sub-cycle: halo update of q;
   - `fv_tp_2d`: `copy_corners`, `yppm` → inner update `q_i`, `xppm` (and the transposed order), flux average;
   - then the flux-form update with `dp2`.
   - Scheme `hord_tr = 12` (PPM with the Huynh-type `pmp/lac` constraint, `lim_fac = 1`) and `nord_tr = 0`
     (no `deln_flux` damping). Confirmed from the run log on 2026-10-16.
3. **Vertical remap** to the hybrid levels: `map1_q2` → `scalar_profile`, `kord_tr = 8`. `fill = F` in
   our runs (checked 2026-10-16), so `fillz` is never called and needs no adjoint.
4. **Global mass fixer** `calcScalingFactorTrop` (`rescale_trop_only = .true.`):
   q ← q · s(q), with s a ratio of global sums. It is nonlinear in q but has a simple rank-1 adjoint.

Around this, in `AdvCore_GridCompMod.F90` / `GCHPctmEnv`: the dry/wet conversions and the pressure
scaling. Their adjoints already exist from the current work. The λ/(air mass) weighting would be
**dropped** in the new mode, because the discrete adjoint acts on λ directly.

Only q (CO2) is an active variable. No wind or pressure adjoint is needed, since the controls are surface
fluxes and, optionally, the initial CO2.

## 3. Existing code we can reuse

**JCSDA/fv3-jedi-linearmodel** (GMAO + JCSDA; files taken from master, last commit d47cea9 of 2022-07-06; GPL v2 headers) contains
Tapenade-generated TL/AD code for exactly these routines (checked 2026-10-16):

| file | routines |
|---|---|
| `tp_core_adm.F90` (9.2 k lines) | `FV_TP_2D_ADM`, `XPPM_ADM`, `YPPM_ADM`, `COPY_CORNERS_ADM`, `PERT_PPM_ADM`, `DELN_FLUX_ADM`, FWD/BWD split |
| `fv_tracer2d_adm.F90` (3.5 k) | `TRACER_2D_FWD/BWD`, `TRACER_2D_1L_FWD/BWD` |
| `fv_mapz_adm.F90` (20 k) | `MAP1_Q2_ADM`, `SCALAR_PROFILE_ADM`, `CS_PROFILE_ADM`, `MAPN_TRACER_ADM`, ... |
| `fv_mp_adm.F90` (3.1 k) | adjoint halo updates `mpp_update_domain2d_*_adm`, `mpp_get_boundary_*_adm`, `fill_corners_*_adm`, group-halo adjoints |

Caveats:
- **Version.** These files are from an older FV3 (no `lim_fac` argument). The `XPPM_ADM` does include
  the `pmp_2/lac_2` constraint branch that `hord 12` uses. A line-by-line diff against our 14.8.0
  `tp_core.F90` / `fv_mapz.F90` is the first task. Version differences would be fixed by re-running
  Tapenade or by hand.
- **The linearization is not exact by default.** JEDI usually runs the TL/AD with linear "pert" schemes
  (`hord_tr_pert = 2`, `split_hord`), not the nonlinear forward scheme. The `_FWD/_BWD` routines do
  support the full scheme with a trajectory. We want that mode (option A below).
- **License.** GPL v2. This is fine for research use. It must be checked before anything is pushed into
  the GCHP/GEOS-Chem forks, which use the GEOS-Chem/MIT-style licenses.

**Tapenade** (INRIA, free for research, Java) can regenerate any routine from our exact 14.8.0 source.
Java 1.8 is on the node. Current Tapenade versions may need a newer JDK, which can be installed in user
space. This is the fallback for routines whose JEDI version is too old.

**TAF** (FastOpt; the user has a license on the laptop, added 2026-10-16) is the preferred generator.
- It produces TL and AD code from our exact 14.8.0 source, which removes the version gap and the GPL
  question: the generated code derives from our own LGPL FV3 source, not from JEDI's.
- It is the tool used for MITgcm/ECCO, so you already know its directives (store/checkpointing,
  `!$taf` flow directives) and its handling of MPI.
- The fv3-jedi-lm files then serve as a reference to compare against, not as the code we ship.
- To check on day 1: whether the license runs locally or through FastOpt's server (the FV3 source is
  open, so sending it is harmless), and how TAF handles FV3's pointers into `gridstruct`. If needed, a
  thin wrapper passes plain arrays.

**Precision.** FV3 in our build is double precision (`FV_PRECISION R8`, `-r8`). The transpose tests must
run in r8 too, so the dot-product identity holds to ~1e-13.

## 4. Options

**A. Exact discrete adjoint of the operational scheme (recommended goal)**
- What it is: the adjoint of `hord 12` / `kord 8` with the limiter branches taken from the forward
  trajectory; the adjoint mass fixer; sub-cycling replayed.
- Exactness: adjoint = transpose of the Jacobian of the forward (to round-off). FD matches except where
  limiter branches switch.
- Needs: the forward q at the start of each advection step. It is inside the step, recomputed (Tapenade
  FWD sweep).

**B. Exact transpose of a linear scheme (stepping stone)**
- What it is: the TL/AD of a linear PPM (`hord 2`/`333`/`-5`) with the forward winds.
- Exactness: exact transpose of a linear operator. Versus the nonlinear forward: limiter errors remain,
  edge/corner asymmetry is removed.
- Needs: no trajectory, so it is simpler.

B is about 60 % of the work of A: the same kernels, halos and integration, without the trajectory. It
can be validated directly with the linear-scheme dot-product test that already exists
(`tests_final/P2_advection/run_P2lin.pbs`, `dotprod_tools.py`). Recommendation: build B first as the
test harness, then switch the kernels to the trajectory-based mode for A.

## 5. Work breakdown

| # | task | details | Claude-assisted (my work + queue) | human alone (experienced adjoint developer) |
|---|---|---|---|---|
| 0 | Audit and design | diff the fv3-jedi-lm vs 14.8.0 routines (`lim_fac`, `hord 12` branch, `nord_tr`, `z_tracer`, `copy_corners`); confirm which checkpointed SPC_ field equals q before advection; define the new `ADV_ADJ_METHOD: discrete` path and unit conventions (λ per kg/kg) | 1–2 days | 1–2 weeks |
| 1 | Kernel adjoints | `xppm/yppm/fv_tp_2d/copy_corners` (+`lim_fac`), `tracer_2d` (sub-cycling, dp2 update), `map1_q2/scalar_profile` (kord 8), mass-fixer adjoint (MPI allreduce) (`fillz` is off, `fill = F`); port from fv3-jedi-lm or regenerate with Tapenade | 3–5 days | 4–6 weeks |
| 2 | Halo / MPI adjoints | adjoint of the q group halo update (accumulate halo values back to the owners, including across face edges with rotation) and of `copy_corners`; from `fv_mp_adm.F90`, adapted to the 14.8.0 FMS/MAPL build | 2–3 days | 2–4 weeks |
| 3 | Standalone unit tests | small driver on C24: TL vs FD, dot-product ⟨TLx,y⟩ = ⟨x,ADy⟩ to ~1e-12 for each kernel and for the full step, on 1, 6, 24 and 48 ranks (decomposition independence), with point seeds at edges/corners | 2–4 days | 2–3 weeks |
| 4 | Integration in AdvCore | new method switch; trajectory from the checkpoints (q before advection, DryPLE0/1, MFX/MFY/CX/CY as now); remove the airmass weighting in this mode; correct order of the backward step relative to the GEOS-Chem adjoint; reduced-checkpoint compatibility | 3–5 days | 3–4 weeks |
| 5 | Validation at C24 (re-use tests_final) | P2a, P2c, P2d/f, T2, P2e, column regression P1, OSSE FD (12 cells), linear-scheme dot test; target: blocks and edges ≲0.5 %, OSSE FD slope 1.000 ± 0.002 | ~1 week wall clock (queue), 1–2 days of my work | 2–3 weeks |
| 6 | C48 check and 4D-Var | repeat P2a/T2 at C48; one E1.3 4D-Var with the new adjoint vs `reverse_wind_pure` (convergence, flux recovery) | ~1 week wall clock (4D-Var run time dominates) | 1–2 weeks |
| 7 | Performance and documentation | timing (expect the advection adjoint ~3–5× the forward advection plus one recompute; advection is a small share of the I/O-bound 4D-Var, so maybe +5–15 % per iteration); memory of the Tapenade stack per rank; notes and paper text | 1–2 days | 1–2 weeks |

### Totals

| | elapsed time | your effort |
|---|---|---|
| **Claude-assisted** | **~4–6 weeks** calendar; about half of it is queue and run time (validation suites, 4D-Var) | **~5–8 working days** spread over that time: design decisions at task 0 (scope A vs B first, license), reviewing the kernel diffs and test results at the end of tasks 3, 5 and 6, approving commits/pushes, and the paper text |
| **Human alone** | experienced FV3/adjoint developer: **~4–6 months full time**; someone new to Tapenade/FV3 internals: 8–12 months | — |

Option B only (stopping before the trajectory): Claude-assisted ~2.5–4 weeks; your effort ~3–5 days.

The estimate assumes the fv3-jedi-lm kernels port with modest changes. If most routines must be
regenerated with Tapenade (larger version gap than expected), add ~1 week (Claude-assisted) or ~1 month
(human).

## 5a. Development kit (C24), ready 2026-10-16

`/home5/ksuselj1/nobackup/adv_adjoint_devkit_C24.tar.gz` (unpacked copy in
`/home5/ksuselj1/nobackup/adv_adjoint_devkit_C24/`; start with its `README.md`). It contains:
- the 14.8.0 FV3 advection sources and GMAO's standalone `tp-core-driver` / `mapz-driver`;
- the fv3-jedi-lm TL/AD reference files;
- per-step C24 advection inputs (MFX/MFY/CX/CY, pressure edges, CO2) for 1 day in three cases (8 m/s,
  tilted rotation, MERRA-2), with initial CO2 and namelists;
- the C24 baselines of the current adjoint (P2a, P2c, P2d, P2e, T2, P2a_lin dot test) with their tools;
- this document and `ADJOINT_PORT_NOTES.md`.

## 6. Risks

- **Limiter non-differentiability.** The exact adjoint is a derivative of a piecewise-smooth map. Small
  FD mismatches at limiter switches are expected and are not bugs. L-BFGS copes with this in practice.
  This would be documented with an ε-scan in the FD checks.
- **Halo adjoint across face edges.** This is the likeliest place for subtle bugs, because of vector
  rotation and corner fills. Task 3 tests it on several decompositions before anything goes into GCHP.
- **Trajectory consistency.** The adjoint must see exactly the forward q, fluxes and `nsplt` of each step.
  Bit-identical reduced checkpoints already exist, but the forward q before advection may not be saved
  at the right point in the step. Task 0 settles this; at worst one extra field per step is checkpointed.
- **Build integration.** The JEDI files depend on their own `fv_arrays_tlmadm` / `fv_mp_adm` modules.
  They will be kept as a separate small library inside FVdycoreCubed, compiled only with `-DADJOINT`.
- **License.** Mentioned in section 3; decide before pushing to the public forks.

## 7. Suggested decision points for you

1. A directly, or B first (recommended: B first; the harness is reused for A).
2. (Updated: TAF is available, so the recommendation is TAF on our own source, with fv3-jedi-lm as the reference.) Use the fv3-jedi-lm code (faster, GPL v2) or regenerate everything with Tapenade/TAF (cleaner license and
   version match, about +1 week).
3. When to start: after the C48 tests/OSSE (the current plan) or in parallel. Most of tasks 0–3 need no
   queue time and could run alongside the C48 jobs.
