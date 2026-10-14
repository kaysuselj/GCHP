# GCHP CO2 adjoint: port 14.5.3 -> 14.8.0 (first, code-only pass)

Written 2026-10-14 (local machine). Source of the port: tag `adjoint-14.5.3-port` on
`kaysuselj/GCHP` (06db738), i.e. the code of `gchp.ref_20261014`. Every change below is
`git diff <upstream 14.5.3 base> <adjoint head>` re-applied by hand / 3-way onto the 14.8.0
release commits; no 14.5.3 branch was merged or rebased.

**Nothing has been compiled against ESMF/MAPL.** The local machine has no ESMF or
netcdf-fortran (14.8.0 needs ESMF >= 8.6.1). Only the Tapenade C runtime (`adStack.c`) was
compiled (clang): it builds and exports all PUSH/POP symbols used by the convection adjoint.
Expect compile errors on Pleiades; they are not hidden bugs, they were never tried.

## 0. Branches

All branches are `adjoint_port_14.8.0`, created from the commit that upstream GCHP 14.8.0
(`e449fcc`) pins for that submodule.

| Component | Fork | 14.8.0 base | Port commits |
|---|---|---|---|
| GCHP (top level) | kaysuselj/GCHP | `e449fcc` (tag 14.8.0) | `e01e3d6` F90 changes; final commit: submodule pointers, `.gitmodules`, this file |
| geos-chem | kaysuselj/geos-chem_adjoint | `a4551f944` | `409904037` |
| FVdycoreCubed_GridComp | kaysuselj/FVdycoreCubed_GridComp | `9d86304` | `f354719` (AdvCore, CMake), `f1d7573` (fvdycore pointer, `.gitmodules`) |
| fvdycore | kaysuselj/GFDL_atmos_cubed_sphere | `643f3c44c` (same as in 14.5.3) | `ca57320` |
| MAPL | kaysuselj/MAPL | `e87909dbc` | `6407a161d` |
| HEMCO | **upstream geoschem/HEMCO 3.13.0, not a fork** | `870e880` | none (see 4.6) |

`.gitmodules` (top level): MAPL, FVdycoreCubed_GridComp, geos-chem point to the forks, branch
`adjoint_port_14.8.0`. `FVdycoreCubed_GridComp/.gitmodules`: fvdycore points to the fork,
branch `adjoint_port_14.8.0`. HEMCO and all other submodules stay upstream.

Not pushed yet (as of writing). Clone after pushing:
`git clone --recursive -b adjoint_port_14.8.0 https://github.com/kaysuselj/GCHP.git`

Build: `-DADJOINT=yes -DREVERSE_OPERATORS=yes` as before. The 14.5.3 top-level
`CMakeLists.txt` hack with hard-coded `/nobackup/ksuselj1/esmf-ESMF_8_2_0` paths was **dropped**
(ESMF 8.2.0 cannot build 14.8.0; use a new ESMF >= 8.6.1 via the normal `ESMF_ROOT`/env). CMake
>= 3.24. The ADJOINT/REVERSE_OPERATORS defines reach GCHP_GridComp and GCHPctmEnv through the
geos-chem build-properties target exactly as in 14.5.3, and FVdycoreCubed through its own CMake
hunk (ported).

## 1. Status per change (handoff section 4)

Legend: **U** applied unchanged, **A** adapted (how), **R** re-derived / needs re-derivation,
**D** dropped (why).

### 1.1 geos-chem

| File / item | Status | Notes |
|---|---|---|
| `GeosCore/Tapenade/` + `GeosCore/CMakeLists.txt` | U (one file D) | 42 of 43 files copied. **D:** `Tapenade/testpushpop` (compiled test executable, binary, cannot be patched and should not be in source). `tapenade_runtime` builds only `adStack.c`, as before. |
| `Interfaces/GCHP/gchp_chunk_mod.F90` | A | Merged cleanly except the `USE HCO_Interface_GC_Mod` line (upstream added `Set_DryDepVel_Diagnostics`; both kept). **AIRQNT adjoint guard changed**, see 2.1. Upstream additions in the forward step (`Set_DryDepVel_Diagnostics` before `DO_MIXING`, `Set_SpcConc_Diags_VVDry`, `PSC2 = PS2`) run in the forward pass of the adjoint step as in a forward run; the adjoint sweep (turbulence, `SfcAdj_PreConv`, convection, surface flux, AIRQNT) is unchanged. |
| `Interfaces/GCHP/Adjoint_Utils_Mod.F90` (new) | U | Dependencies checked against 14.8.0: `Convert_Spc_Units` (now in `GeosUtil/unitconv_mod.F90`, same interface), `KG_SPECIES_PER_KG_DRY_AIR`, `MAPL_CommsAllReduceMin/Sum`, `MAPL_CommsBcast` (MAPL 2.69 overloads exist). |
| `Interfaces/GCHP/Chem_GridCompMod.F90` | A | Merged cleanly; three adaptations: (a) `DELP_DRY` restore from the checkpoint now on **every** adjoint call incl. the first (2.2); (b) `ADV_AIRDEN` export spec and its fill now inside `#ifdef ADJOINT` (the AdvCore import was ADJOINT-only already; unguarded, a non-ADJOINT build would fail the connectivity); (c) the debug `print *,'Kay before chunk run',Species(3)%Conc(6,5,1)` was unguarded (printed in every build, out of bounds on small domains), now `#ifdef ADJOINT`. `IsAdjTime` alarm, checkpoint read/write, SpeciesAdj refresh, `Adjoint_RefreshExports`, forcing at `t+dt`: unchanged. |
| `GeosCore/convection_mod.F90` | A (no re-derivation needed for CO2) | See 2.3. Adjoint block (forward replay with PUSH/POP + reverse sweep) moved into the renamed `DO_RAS_CLOUD_CONVECTION`, after upstream's new washout initialisation. New: stop with an error if an adjoint run selects Grell-Freitas convection (`convection_mod.F90:141`). |
| `GeosCore/pbl_mix_mod.F90` | A | TurbDay adjoint unchanged numerically. Moved the adjoint setup block (it sat between `!$OMP PARALLEL DO` and its `DO`, invalid OpenMP) above the directive; added `ADJS` to `PRIVATE`; fixed `#ifdef adjoint` -> `#ifdef ADJOINT` (pointer nullify was never compiled). Serial behaviour identical. Upstream did not change this file. |
| `GeosCore/mixing_mod.F90` | U | The 14.5.3 diff only removes debug prints (no reordering in this diff, contrary to the handoff wording). Upstream changed one line. The upstream ADJOINT branches in `DO_TEND` (negative TS, drydep `FRAC` on SpeciesAdj) are unchanged since 14.5.3. |
| `GeosCore/hco_interface_gc_mod.F90` | A | `Compute_Sflx_For_Adjoint`, `Get_Adjoint_Hemco_SurfaceFlux` added unchanged (self-contained, `#ifdef ADJOINT`); PUBLIC line conflict resolved (kept upstream `Set_DryDepVel_Diagnostics`). HEMCO APIs checked in 3.13.0: `HcoDiagn_AutoUpdate`, `Diagn_Create(... AutoFill, COL, OkIfExist)`, `GetExtNr`, `HcoState%Diagn%HcoDiagnIDDefault`, `HcoState%Config%ExtList`, `HCO_GC_GetDiagn` (2D), `GetHcoValEmis`, `LoadHcoValEmis`, `InquireHco`. The new EFLX/DFLX debug WRITE inside `Compute_Sflx_For_Vdiff` was unguarded (forward log change), now `#ifdef ADJOINT`. |
| `GeosCore/hco_utilities_gc_mod.F90` | U | Only null-pointer guards in two `HCO_GC_GetDiagn` routines. The "scale-factor options 6 and 7" of the handoff are not in this file: they are `FD_STEP: 6/7` (perturbation x1.01 / x0.99) in `Setup_Adjoint_ForwardPert` (`Adjoint_Utils_Mod.F90`), ported unchanged. |
| `Headers/state_chm_mod.F90`, `state_diag_mod.F90`, `input_opt_mod.F90`, `GeosCore/diagnostics_mod.F90`, `GeosCore/input_mod.F90`, `Interfaces/GCHP/Includes_Before_Run.H`, `Interfaces/GCHP/CMakeLists.txt` | U | Clean 3-way merge. `CF_IMIN...` fields commented out; checked that no live 14.8.0 code uses them. `input_mod.F90` comments out three log WRITEs (`P_ICE_SUPERSAT`, `GAMMA_HO2`, `AUTOREDUCE_THRESHOLD`) in all builds: log only, kept 1:1. |
| `CHANGELOG.md` | D | Documentation only. |

### 1.2 GCHP top level

| File / item | Status | Notes |
|---|---|---|
| `GCHPctmEnv_GridCompMod.F90` (a) first-step advection fix | U (intent) | Upstream 14.8.0 still skips `fv_computeMassFluxes` on the first call in every ADJOINT build. Fixed as before: skip only when `dt < 0` (`GCHPctmEnv_GridCompMod.F90:1508`). Also see 2.4: AdvCore had the same skip. |
| `GCHPctmEnv_GridCompMod.F90` (b) `ADV_ADJ_MASSFLUX: forward` | U | Reverse run: Courant numbers of the reversed wind, mass fluxes `= -fv_computeMassFluxes(..., -dt)`. Module flag now inside `#ifdef ADJOINT` (MAPL2 section). The `dt before prepare_massflux_exports` print was unguarded in 14.5.3, now `#ifdef ADJOINT`. Only the wind path is covered; with `IMPORT_MASS_FLUX_FROM_EXTDATA: .true.` the reverse run uses ExtData fluxes as before (untested, as in 14.5.3). |
| `GCHP_GridCompMod.F90` | A | `Adjoint_RefreshExports` after the reversed advection: unchanged. `ADV_AIRDEN` connectivity now `#ifdef ADJOINT` (see Chem_GridCompMod). Placed after upstream's new `DELPDRY` CHEM->ADV connectivity. |
| `CMakeLists.txt` | D | The 14.5.3 diff here was not the ADJOINT options but hard-coded Pleiades ESMF 8.2.0 paths. |

### 1.3 FVdycoreCubed_GridComp

| File / item | Status | Notes |
|---|---|---|
| `AdvCore_GridCompMod.F90` | A | 3-way merge, 6 conflicts resolved by hand. Ported unchanged in behaviour: `ADV_ADJ_METHOD` guard (`reverse_wind_pure` only), `ADV_ADJ_WEIGHT` (`dpadv` default, `airmass`, `airden`), `ADV_ADJ_RESCALE`, `ADV_AIRDEN` import, `PLE0/PLE1` and `DryPLE0/DryPLE1` swapped in reverse runs, `advOn` (no advection at the clock stop time in reverse runs; `AdvCore_Advection>0` otherwise), dpadv division before / multiplication with `PLEAdv` + `(dp(t)-dp(t+dt))*adjMuIn` after, `fill` off for adjoint runs during `offline_tracer_advection`, post-advection block only if `advOn`, global-sum rescale. Adapted: upstream now clears `firstRun` at the end of `Run` (was: before the advection call); the advection guard is `advOn` in both ADJOINT and non-ADJOINT builds, nothing between the advection call and the end of `Run` reads `firstRun`, so behaviour is as in 14.5.3. Upstream 14.8.0 still had `AdvCore_Advection>0 .and. .not. firstRun` in ADJOINT builds, i.e. the forward first-step skip existed in AdvCore too; it is fixed by `advOn`. **New upstream DELPDRY start-up scaling skipped in adjoint runs** (2.5). Dropped: unused `IMPORT_MASS_FLUX_FROM_EXTDATA` resource read (only used by the removed `reverse_wind` variant). |
| `FV_StateMod.F90` | D | The 14.5.3 diff is **not** "hooks / direction of MFX, MFY, PLE0, PLE1" as the handoff says: it only adds `Log_AdjCO2_Dynamics_Sum` in `FV_Run` (full dynamics, not used by GCHP's offline advection) whose output is commented out (it does an allreduce and nothing else) plus unused locals. No effect on results. |
| `CMakeLists.txt` | U | ADJOINT/REVERSE_OPERATORS compile definitions (upstream still does not define them for this target). |

### 1.4 fvdycore

| Item | Status | Notes |
|---|---|---|
| `model/fv_tracer2d.F90` halo zeroing of `dpL`, `q3` | U | fvdycore is at the same commit `643f3c44c` in 14.5.3 and 14.8.0; the hunk applied verbatim. |

### 1.5 MAPL

| Item | Status | Notes |
|---|---|---|
| `MAPL_CapGridComp.F90` `step_reverse`: ExtData after the clock step-back | U | `step_reverse` in the 14.8.0 MAPL (2.69 line) is identical to the 2.26 one; hunk applied with offset. The reverse-time loop, `REVERSE_TIME:` handling and History `sticky=.false.` changes are all present upstream. |

### 1.6 HEMCO

No port; upstream 3.13.0 kept. The only adjoint-era source change in 14.5.3
(`CALL HCO_MSG(SEP1=' ', ...)` in `hco_calc_mod.F90`) no longer occurs in 3.13.0.

## 2. Things to look at first (unsure / changed semantics)

### 2.1 AIRQNT mixing-ratio update (14.8.0 change)
14.8.0 `AirQnt` defaults `update_mixing_ratio` to false, and updates the mixing ratios only if
the flag is true **and** `|sum(DP_DRY_PREV) - sum(DELP_DRY)| < 1e-14` (`calc_met_mod.F90`). The
chunk still passes `scaleMR` (true from the 2nd step if `LTRAN` is off). The 14.5.3 adjoint was
guarded by `scaleMR` only. Now: `AirQntUpdatedMR` = the flag passed to that step's AirQnt call
AND the same global-sum test, evaluated right after the call (`gchp_chunk_mod.F90:1133`), guards
the adjoint (`gchp_chunk_mod.F90:1992`). With `LTRAN` on (normal case) both are false, as before.
The 1e-14 absolute test on sums of hPa values looks odd upstream (true essentially only when the
pressures are identical); whatever it does, the adjoint now mirrors it.

### 2.2 DELP_DRY restore in adjoint mode
14.5.3 upstream restored `State_Met%DELP_DRY` (and BXHEIGHT, TropLev) from the internal state on
the first call; the adjoint added the restore on every later adjoint call (`.NOT. FIRST`).
14.8.0 removed the upstream first-call block (restart DELP_DRY now goes to FV3 as `DELPDRY`).
The adjoint restore is now unconditional in adjoint runs (`Chem_GridCompMod.F90:3072`), which
gives the adjoint the same DELP_DRY as in 14.5.3 also on the first call. BXHEIGHT/TropLev are no
longer restored on the first call (upstream change, affects both modes, not adjoint-specific).

### 2.3 Convection (re-derivation scoping)
Compared upstream RAS `DO_CLOUD_CONVECTION` (14.5.3) vs `DO_RAS_CLOUD_CONVECTION` (14.8.0):
- **Unchanged for an insoluble species:** cloud base search (`DQRCU>0`), below-cloud mixing
  (`QB`, `QC`, `MB`), in-cloud mass balance (`ENTRN`, `CMOUT`, `QC_PRES/SCAV`, `T1..T4`, `DELQ`,
  negativity clip), internal time-step loop, `BMASS`.
- Changed but irrelevant for CO2 (adjoint assumes F=0, no washout, no Hg): washout restricted to
  `NW>0` and to `LUO_WETDEP`; `PFICU/PFLCU` slices `2:NZ+1` -> `1:NZ` (only `PDOWN`/`ALPHA`);
  GF precip reconstruction removed from RAS (`DQRCU`, `REEVAPCN` now direct met pointers; same
  values for RAS met); `QC = Q(CLDBASE)` set before the branch (equivalent); T0..DELQ zeroed.
- Therefore the existing Tapenade-style adjoint is still the exact transpose for CO2; **no
  re-derivation was done or needed**. Its forward replay uses the local `PDOWN` built from the
  new `1:NZ` slices, consistent with the new forward.
- GF (`Grell_Freitas_Convection`, GEOS-IT/FP) has no adjoint; adjoint runs stop with an error.
- Pre-existing, not changed (port 1:1): uses `NA = NFD` both as species and advected index
  (`Map_Advect(NA)`, `F(k,na)`, `SpeciesAdj(...,NA)`): correct only if CO2 is species 1 and
  advected species 1 (carbon CO2-only); PUSH/POP use a global stack and `DO_CONVECTION` has an
  OpenMP loop over columns: requires OpenMP off (as in 14.5.3); `PUSHREAL4ARRAY(x, fp/4)` hack.

### 2.4 First-step advection
Fixed in both places where upstream skips step 1 in ADJOINT builds: GCHPctmEnv mass fluxes and
the AdvCore `.not. firstRun` guard. Verify on Pleiades that a forward ADJOINT build now matches a
non-ADJOINT build bit for bit from step 1 (best single regression test of this port).

### 2.5 New first-step DELPDRY scaling in AdvCore (14.8.0) - adjoint missing
On the first AdvCore call with `DELPDRY` from the restart (`sum > 0`), 14.8.0 scales all advected
tracers by restart dp / met dp (`scale_tracers_by_pressure_ratio`). In a reverse run the first
call is at the end of the window, so applying it there would be wrong; it is **skipped in adjoint
runs** (`AdvCore_GridCompMod.F90:1258`). The transpose of the forward scaling (multiply the
adjoint by the same per-cell ratio after the last advection adjoint, clock = window start) is
**not implemented (TODO)**. Until then: either use restarts without `DELP_DRY` / with zero
`DELP_DRY` for forward/FD runs (then the forward does not scale and the adjoint is consistent),
or check that the ratio is 1 (restart written by the same met). FD tests of the initial-CO2
adjoint will show it as an error of size |ratio - 1| otherwise.

### 2.6 Met timing (14.8.0 change)
`State_Met%T`, `SPHU` are now `TMPU2`, `SPHU2` (post-advection, was the average of 1 and 2), and
`PSC2_WET/DRY` start from `PS2` (was `PS1`) before `SET_FLOATING_PRESSURES`. The adjoint itself
does not use T/SPHU; it uses them only through the forward replay of the step (AIRDEN ->
`ADV_AIRDEN` for `ADV_ADJ_WEIGHT: airden`, `AD`/`DELP_DRY` in AIRQNT, convection `BMASS` via
`DELP_DRY`, PBL `AD`). The forward pass of each adjoint step reads the same ExtData fields as the
forward step (MAPL fix), so these should still agree; the step-by-step forward-vs-adjoint met
comparison (acceptance test 3) must be repeated with the 14.8.0 field list, including `PS2`,
`SPHU2`, `TMPU2`.

### 2.7 Other remarks
- `Set_DryDepVel_Diagnostics` (new, before `DO_MIXING`) changes drydep velocities via HEMCO
  sea-air exchange; CO2 has no drydep in the carbon sim, so no adjoint impact expected.
- AdvCore sets `FV_Atm(1)%flagstruct%fill = .true.` after the adjoint advection call
  unconditionally (as in 14.5.3), not the previous value; fine if `fill` is true in the namelist.
- `ADVCORE_ADJ_DEBUG` debug blocks and the `ADJ_CO2_SUM` allreduce prints in
  Chem_GridCompMod were kept 1:1 (they do allreduces every step in adjoint runs).
- History `SpeciesAdj` refresh (`Adjoint_RefreshExports`) kept; end-time adjoint record still
  zero; vertical-order and units pitfalls of handoff section 7 unchanged.

## 3. OSSE / satellite-track sampling (not ported)

MAPL 2.69 History has no `track_file` sampler (present in 14.5.3's MAPL History). It has
`sampler_type: trajectory` (`gridcomps/History/Sampler/MAPL_TrajectoryMod*.F90`) and a station
sampler. Needed for the OSSE: (1) convert OCO-2 track files to the trajectory input format
(level dimension first, see MAPL docs); (2) check that the trajectory sampler works with a
reversed clock / backward collections (it is time-windowed; never tested in reverse); (3) adapt
the HISTORY.rc collections and the post-processing that reads the sampled output. The adjoint
forcing loader (`Load_OCO2_Adjoint_Forcing`) does not use History and is unaffected.

## 4. Suggested test order on Pleiades

1. Build (ADJOINT, REVERSE_OPERATORS, new ESMF). Non-ADJOINT build too (guards above).
2. Forward ADJOINT build vs non-ADJOINT build, bit-identical from step 1 (2.4).
3. Forward vs adjoint met/emissions step by step (2.6), forcing loader 6 vs 24 cores.
4. P1 column tests (advection off): initial-CO2 adjoint vs FD (2.1, 2.5), surface-flux vs FD
   with PBL + convection (2.3).
5. T6 / P2 advection tests, then P3.
