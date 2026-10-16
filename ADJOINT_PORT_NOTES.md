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

## 5. Pleiades results

Session on Pleiades from 2026-10-14 (layout: `gchp_14.8.0_adj/code` = this clone, `gchp_14.8.0_adj/tests_final`).

### 5.1 Build environment
- No ESMF module on Pleiades: ESMF 8.6.1 built from source (`/nobackup/ksuselj1/esmf-ESMF_8_6_1_intel2023`,
  PIO/LAPACK internal, netCDF via nc-config). The global `ESMF_SITE=NAS` of the Pleiades environment needs a
  `build_config/Linux.intel.NAS` directory (copy of `Linux.intel.default`) for `make install`.
- Compiler: Intel 2023.2.1 (`comp-intel/2023.2.1`, ifort 2021.10), HPE MPT 2.30, netCDF-C 4.4.1.1 /
  netCDF-Fortran 4.4.4 (`netcdf/4.4.1.1_mpt`), CMake 3.26.5. Environment: `gchp_14.8.0_adj/gchp_14.8.0_env.sh`
  (load the modules one by one: a single `module load` line with all of them silently loads nothing).
- MAPL 2.69 needs udunits2: `/nasa/gdl/udunits_2.2.28`, passed as
  `-Dudunits_LIBRARY=.../lib/libudunits2.so -Dudunits_INCLUDE_DIR=.../include -Dudunits_XML_PATH=.../share/udunits/udunits2.xml`.
- Configure: `cmake .. -DMECH=carbon -DADJOINT=yes -DREVERSE_OPERATORS=yes <udunits options>` in `build_adj`.

### 5.2 Compile fixes
| # | File | Fix | Why |
|---|---|---|---|
| 1 | MAPL `gridcomps/History/Sampler/MAPL_MaskMod.F90` (upstream MAPL code, not adjoint) | removed the two `use pflogger, only: Logger, logging` inside the separate-module-procedure interface bodies (`MaskSampler_from_config`, `create_Geosat_grid_find_mask`); `Logger`/`logging` come from the host module's own `use pflogger` | internal compiler error in ifort 19.1 (2020.4) and 2021.10 (2023.2.1), also at -O0; ifx 2023.2 internal error; ifx 2025.3 rejects it (error #6405, same named entity LOGGER). No behaviour change. |
| 2 | geos-chem `GeosCore/convection_mod.F90` (`DO_RAS_CLOUD_CONVECTION`) | the adjoint `ELSE` branch of `IF (is_forward_run)` (from `NA=Input_Opt%NFD` to the end of the reverse sweep) is now inside `#ifdef ADJOINT ... #endif` | the plain (non-ADJOINT) build failed: the branch uses `QB0`, `branch`, `ad_count`, ... (declared only under `#ifdef ADJOINT`), `Input_Opt%NFD` and `State_Chm%SpeciesAdj`. ADJOINT build: code unchanged. |
| 3 | geos-chem `Interfaces/GCHP/Adjoint_Utils_Mod.F90` | whole module inside `#ifdef ADJOINT` (empty object in plain builds) | plain build: uses ADJOINT-only `Input_Opt` fields (`NFD`, `IFD`, `IS_FD_*`, ...); its only user (`gchp_chunk_mod.F90`) already has the `USE` under `#ifdef ADJOINT`. |
| 4 | geos-chem `Interfaces/GCHP/gchp_chunk_mod.F90` (`GCHP_PRINT_MET`, debug print) | the two WRITEs of `SurfaceFluxAdj`/`SpeciesAdj` inside `#ifdef ADJOINT` | plain build: ADJOINT-only `State_Chm` fields. |
| 5 | MAPL `gridcomps/History/MAPL_HistoryGridComp.F90` (end_alarm for collections without end_date) | backwards collections: `end_alarm` created disabled (`enabled=.false., sticky=.false.`), as in the GCHP MAPL 2.26 | MAPL 2.55+ creates it at CurrTime (= window end in a reverse run); it rings at the last reverse step ("Not Writing alarm ... because end_alarm is ringing") so the adjoint record at the window START (the initial-condition sensitivity) was never written. The MAPL source itself flags this ("Might need to be adapted for GCHP adjoint"). Forward collections unchanged (T0 stays valid). Binary `build_adj/bin/gchp.adj_20261015`. |

### 5.3 Builds
- ADJOINT build (`build_adj`): compiles and links with fix 1 only; adjoint code is in the binary (ADV_ADJ_* keys,
  Load_OCO2, Tapenade PUSH/POP symbols). Frozen as `build_adj/bin/gchp.adj_20261014` (built before fix 2, which
  does not change ADJOINT builds).
- Plain build (`build_plain`): needed fix 2.

- Plain build (`build_plain`) compiles with fixes 1-4: `build_plain/bin/gchp.plain_20261014`. The ADJOINT build
  was recompiled from the same source (fixes 2-4 do not change ADJOINT code): `build_adj/bin/gchp.adj_20261014`.

### 5.4 Run directory (tests_final/rundir_carbon_template)
- `createRunDir.sh`: carbon, all species (as 14.5.3: CO2 = species 3 = advected 3, so the `NA = NFD` caveat of 2.3 does
  not apply), MERRA-2. Met: 14.8.0 reads the same 73 MERRA-2 fields as 14.5.3 (only PS2/SPHU2/TMPU2 read frequency
  follows the dynamic step); MetDir -> the 4x5 archive `/nobackupp17/lwu5/GCHP/ExtData/GEOS_4x5/MERRA2`, file names
  `.05x0625.` -> `.4x5.` in ExtData.rc. All met files for Jan-Mar 2016 present.
- HEMCO inputs from `GC_DATA_ROOT=/nobackup/ksuselj1/geos_chem_new/ExtData`: 215 of 225 ExtData entries present;
  missing GFEIv3 / GRPI (new CH4 inventories) and GFAS. Switched off GFEIv3, GRPI, GFAS; GFED on (fires). No effect on CO2.
- `tests_final/common/make_rundir.sh`: C24, 24 cores, 1800 s steps, start 2016-01-01, restart link.
- Restarts: `tests_final/common/GEOSChem.Restart.20160101_0000z.c24.nc4` (14.5.3 restart, with DELP_DRY) and
  `...c24.noDELP.nc4` (DELP_DRY = 0, for adjoint/FD tests, notes 2.5).
- Batch jobs: non-login shells have no `module` command; the env file now sources `/usr/share/Modules/init/bash`;
  14.8.0's `setCommonRunSettings.sh` and `setRestartLink.sh` must be sourced without `set -eu`.

### 5.5 Test results
| Test | Setup | Result |
|---|---|---|
| T0 forward ADJOINT build vs plain build | CO2-only carbon run dir (`rundir_co2_template`), C24, 24 cores, 1 day, 1800 s, advection + convection + PBL mixing on, chemistry/deposition/non-local mixing off, restart with DELP_DRY (start-up scaling active in both) | **PASS: bit-identical**, 24 hourly SpeciesConc (instantaneous) and 216 emission fields (`tests_final/T0_adj_vs_plain/T0_report.txt`) |
| P1 column PBL mixing (`pbl`) | `tests_final/P1_columns/run_P1.pbs`, C24 CO2-only, 1 day, 1800 s, advection off, PBL mixing on; fluxes as the 14.5.3 tests/OSSE (template switched 2026-10-15: FF_regrid2, ECCO-Darwin, NBE CARDAMOM, GPP ORCHIDEE, TER CLASS-CTEM; sigma 750 on TER); conc seeds LAYER r=1,10,20,30; FD from 14.8.0 forward runs (`gchp.adj_20261015`) | **PASS**, same as 14.5.3: conc median 1.000000, max abs(x-1) 1.0e-6; flux (r=1,10, N=1324/624 as 14.5.3): median 1.000011 after MW_CO2/MW_air, max abs(x-1) 1.2e-5; FD surface response 4.242e-5 (14.5.3: 4.243e-5) |
| P1 column PBL + convection (`pblconv`) | as above with convection on | **PASS**, same as 14.5.3: conc r=1..30 max abs(x-1) 1.2e-6; flux r=1..30 (N=1324/1014/498/318; 14.5.3 1324/1015/497/318) median 1.000011, max abs(x-1) 1.2e-5; FD surface response 4.179e-5 (14.5.3: 4.180e-5) |
| P1c forcing from file, 6 vs 24 cores / P1d forcing timing | `tests_final/P1_columns/run_P1cd.pbs` on the pblconv setup: OSSE forcing files of 2016-01-01 (copied from E1.3); synthetic level-10 forcing file stamped at the window end vs the LAYER seed r=10 | **PASS, as 14.5.3**: 6 and 24 cores bit-identical (170 obs read, 170 deposited); file seed = built-in seed bit-identical in all 49 records (`pblconv/results_P1cd/P1c_compare.txt`, `P1d_compare.txt`) |
| P2a face edge/corner, 8 m/s zonal wind | `tests_final/P2_advection/run_P2.pbs` (EXP=P2a): advection only, C24, 5 days, idealized met (MetIdealized), restart DELP_DRY=0; receptors, seeds and the 57 perturbed cells of the 14.5.3 `_ref` experiment; FD from 14.8.0 forwards; ADV_ADJ_METHOD reverse_wind_pure, airmass, massflux forward | **Same as 14.5.3.** FD values and J_base identical to 14.5.3 to >= 5 digits (forward advection unchanged for these winds). Own perturbations n=52: median |rel err| 21.8 % (14.5.3: 21.8 %), rms 164.7 % (163.9 %), corr 0.9969 (0.9969), sum ratio 1.0269 (1.0267). Per receptor (median, 14.8.0 / 14.5.3): westofseam 1.8/1.9 %, onseam 22.4/22.3 %, eastofseam 8.3/8.1 %, corner_west 50.4/50.6 %, corner_on 37.1/37.0 %, corner_east 80.8/80.3 %. The face-edge/corner adjoint error of 14.5.3 is unchanged in 14.8.0. Report `P2_advection/P2a/results/report.txt` |
| P2c tilted rotation | `run_P2.pbs` EXP=P2c (airmass) and EXP=P2c_dpadv (ADJ_WEIGHT=dpadv), MetIdealized_TiltedRot, 20 cells of the 14.5.3 `_ref` experiment | FD identical to 14.5.3. inside_f1: median 4.5 % (airmass) / 3.5 % (dpadv; 14.5.3 3.5 %). cross_f2_f1 (face 2 -> 1 near 90E, Tibet margin): airmass median 46.1 %, corr 0.79, sum 0.81; **dpadv median 26.8 %, corr 0.984, sum 0.918 = 14.5.3 exactly** (the 14.5.3 `_ref` run used the binary default dpadv). No port difference; the weighting choice matters at this crossing (dpadv better here, airmass better in the real-wind block tests T2/T5b of 14.5.3). |
| P2d real winds, 5 days (+ P2f) | `tests_final/P2_advection/run_ADV.pbs` EXP=P2d: MERRA-2 4x5, 6 point receptors, 18 patch perturbations of the 14.5.3 `P2d_real_winds_5days_ref_airmass`, airmass | **Identical to 14.5.3** to all printed digits for every receptor (own perturbations: median 25.4 %, corr 0.905, sum 0.848; 14.5.3 the same). P2f: forcing-file seed vs built-in REGIONAL seed bit-identical in all 241 records (as 14.5.3) |
| T2 block receptors, real winds, 5 days | `run_ADV.pbs` EXP=T2 (tadv_tools), receptors/perturbations of 14.5.3 `T2_block_receptor_real_winds_ref` | **Same as 14.5.3**: weighted error tropical block 4.4 % (14.5.3 4.5 %), jet block 5.5 % (5.2 %); point receptors 45.1 % / 33.5 % (44.5 / 34.9 %). Resolved-scale validation holds in 14.8.0 |
| P2e real winds, 1 month, 6 block receptors | `run_ADV.pbs` EXP=P2e (tadv_tools), receptors/15 perturbations of 14.5.3 `P2e_real_winds_1month_blocks_ref`, ADV_ADJ_WEIGHT airmass | Weighted error jet 4.0 %, tropical 5.7 %, edge 7.3 %, corner_N 8.2 %, corner_S 12.9 %, polar 7.3 %; sum ratios 0.92-1.01. The 14.5.3 `_ref` run used the binary default **dpadv** (21-67 %), so it is not a like-for-like reference; with airmass the month-long block validation is new and holds at 4-13 % |
| OSSE full-physics gradient (replaces P3b) | 4dvar_final E1.3_const_sat_148 (truth with 14.8.0, offline satellite sampler, 3 months, transport+convection+PBL): central FD (eps 0.05) of J_obs vs iteration-1 adjoint dJ/dsigma at 12 cells (`fd_check/run_fd148.pbs`) | **PASS, identical to 14.5.3**: slope 0.9940, corr 0.9993, weighted error 1.07 %, cells 0.961-1.009 (14.5.3: 0.9940 / 0.9993 / 1.07 %); iteration-1 gradient vs 14.5.3 corr 1.000000 |

Run notes:
- Both builds abort in finalize with `double free or corruption` AFTER the last step and the restart write
  (MAPL/ESMF/MPT finalize on Pleiades; not adjoint-specific; output and checkpoint are complete). Test scripts accept a
  run that reached its end time and wrote `gcchem_internal_checkpoint`. Not yet investigated.
- One early 3-day attempt hung before the first step (0 steps in 2 h); not reproduced since.
- 14.8.0 run dirs: History frequencies from `setCommonRunSettings.sh` apply only with `AutoUpdate_Diagnostics=ON`.
- Upstream 14.8.0 `setCommonRunSettings.sh` (adjoint branch, ~line 737): `sed "s/^\([\t ]*\)'Adjoint',/\1'Adjoint',/"`
  is a no-op, so an `#'Adjoint',` collection is never re-enabled for `Model_Phase=ADJOINT` (same for `SFEmissions`).
  Tests enable it directly. Should be `s/^\([\t ]*\)#'Adjoint',/\1'Adjoint',/`.
- Reverse (adjoint) runs print no "GCHP Date" progress lines in 14.8.0; they finish with the same finalize
  `double free` after writing the final checkpoint.
- P1: the old `column_fd_tools.py analyze` (mass-ratio regression) crashes on zero-mass columns (SVD not converged); the
  `column_checks_v2.py` checks above are the pass criterion (same as 14.5.3). Before fix 5 the window-start (r=1 IC) record was missing.

- P1 first run (stock 14.8.0 CO2 fluxes, sigma on Net_terrestrial_exch) gave the same ratios but a 9x smaller flux response; results kept in `P1_columns/stockflux_20261015/`.

- P2 run notes: the first adjoint step reports one "Missing adjoint checkpoint file ...<window end>" + MAPL FAIL (status=-1) per rank; also in 14.5.3, harmless (run_P2.pbs allows exactly this).
- 2026-10-15 /nobackupp27 outage (MDT failover, ~10:00-16:00) truncated some files written with sed -i in run dirs to 0 bytes (P2c adjoint GCHP.rc/input.nml, fd_w1 ExtData.rc); dirs rebuilt from the template.



- 14.8.0 stock ExtData reads the (unused) ODIAC weekly factors with F2006-01-%d2 -> bracket failure on day 31 (any run crossing a 31st). Fixed in the CO2 template/OSSE dirs by loading factors 440/441/480 only for FOSSIL_CDIAC/CO2CORR (both off); no flux change.


