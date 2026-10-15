# Handoff: build and test the GCHP 14.8.0 adjoint port on Pleiades

Written 2026-10-14 on the user's local machine for a Claude Code session on Pleiades (NASA NAS).
Author of the adjoint: Kay Suselj (GitHub `kaysuselj`). Read this whole file, then
`ADJOINT_PORT_NOTES.md` (same directory), before touching code.

## 1. Task

The CO2 adjoint (and its 4D-Var hooks) was ported from GCHP 14.5.3 to GCHP 14.8.0 in a
code-only pass. **It has never been compiled against ESMF/MAPL** (no ESMF on the local
machine). Your job:

1. Build the port on Pleiades (ADJOINT and non-ADJOINT), fixing compile errors.
2. Run the test sequence in section 6 and compare with the 14.5.3 reference numbers.
3. Record every fix and every result in `ADJOINT_PORT_NOTES.md` (add a section
   "Pleiades results"), so the next session can continue from it.

Compile errors are expected and are not a sign that the port is wrong in concept. Fix them by
keeping the intent described in `ADJOINT_PORT_NOTES.md`. If a fix changes numerics or forward
behaviour, stop and ask the user first.

## 2. Code

| Component (path in GCHP) | Fork | Branch | 14.8.0 base | Port head (2026-10-14) |
|---|---|---|---|---|
| GCHP (top level) | kaysuselj/GCHP | `adjoint_port_14.8.0` | `e449fcc` (tag 14.8.0) | this commit |
| `src/GCHP_GridComp/GEOSChem_GridComp/geos-chem` | kaysuselj/geos-chem_adjoint | `adjoint_port_14.8.0` | `a4551f944` | `409904037` |
| `src/GCHP_GridComp/FVdycoreCubed_GridComp` | kaysuselj/FVdycoreCubed_GridComp | `adjoint_port_14.8.0` | `9d86304` | `f1d7573` |
| `.../FVdycoreCubed_GridComp/fvdycore` | kaysuselj/GFDL_atmos_cubed_sphere | `adjoint_port_14.8.0` | `643f3c44c` | `ca57320` |
| `src/MAPL` | kaysuselj/MAPL | `adjoint_port_14.8.0` | `e87909dbc` | `6407a161d` |
| HEMCO | upstream geoschem/HEMCO 3.13.0 | (not a fork) | `870e880` | no changes |

Clone:

```
git clone --recursive -b adjoint_port_14.8.0 https://github.com/kaysuselj/GCHP.git gchp_14.8.0_adj
```

`.gitmodules` already points MAPL, FVdycoreCubed_GridComp (and its fvdycore) and geos-chem to
the forks. Check after cloning: `git submodule status --recursive` must show the port heads
above. The complete port of a component is `git diff <14.8.0 base> <port head>` in that repo.

Reference (14.5.3): tag `adjoint-14.5.3-port` on kaysuselj/GCHP, binary `gchp.ref_20261014`
(and `gchp.dpadv_20261014`, the paper binary). Use these for comparison; do not modify them.

## 3. Build

- Needs **ESMF >= 8.6.1** and **CMake >= 3.24**. The existing 14.5.3 build uses ESMF 8.2.0
  (`/nobackup/ksuselj1/esmf-ESMF_8_2_0`), which cannot build 14.8.0; a new ESMF build is
  needed. Ask the user which ESMF/compiler/MPI modules to use if it is not obvious; do not
  overwrite the existing 8.2.0 install.
- The 14.5.3 trick of hard-coding the ESMF paths in the top-level `CMakeLists.txt` was
  removed. Point CMake at ESMF the standard GCHP way (`ESMF_ROOT` / `ESMFMKFILE`
  environment).
- 14.8.0 merged the support libraries (gFTL, pFlogger, yaFyaml, pFUnit, fArgParse) into the
  `src/GFE` submodule; they are built in-tree.
- Configure two build directories from the same source:
  - adjoint: `cmake <src> -DADJOINT=yes -DREVERSE_OPERATORS=yes -DRUNDIR=<rundir>`
  - plain: `cmake <src> -DRUNDIR=<rundir>` (standard GCHP; needed for test T0)
- Carbon simulation, MERRA-2 met, c24, as in 14.5.3. OpenMP must be off (the convection adjoint
  uses the global Tapenade PUSH/POP stack inside the column loop of `DO_CONVECTION`).
- Build the **plain** configuration too: several adjoint additions are now behind
  `#ifdef ADJOINT` and a plain build checks that nothing leaked out.

### Likely compile problems (none of these were seen; they are guesses)
- `gchp_chunk_mod.F90`: calls `Print_Global_Species_Kg(...)` and `GCHP_PRINT_MET(...)`
  (signature extended with `State_Chm`); check argument lists against 14.8.0
  `GeosUtil/unitconv_mod.F90`.
- `convection_mod.F90`: the adjoint block in `DO_RAS_CLOUD_CONVECTION` uses lower-case names
  of locals (`dqrcu`, `pflcu`, `pficu`, `reevapcn`, `pdown`, `bmass`) that are pointers or
  arrays in 14.8.0; Tapenade calls (`PUSHREAL4ARRAY`, `PUSHCONTROL1B/2B`, `PUSHINTEGER4` and
  the POP versions) are external, provided by `GeosCore/Tapenade/adStack.c`. Module-level
  `USE MAPL_CommsMod, only: MAPL_AM_I_ROOT` is unguarded.
- `Chem_GridCompMod.F90`: new locals (`ForceAdjointInitRun`, `AdjLocalSum`, `VM` use in
  `MAPL_CommsAllReduceSum`, `tsChemInt`, `adj_elapsed_s`); check they exist in the 14.8.0
  routine scope.
- `AdvCore_GridCompMod.F90`: many new locals in `Run`; `MAPL_CommsAllReduceSum` needs
  `MAPL_Mod` (already used).
- Link: `tapenade_runtime` is linked `PRIVATE` into `GeosCore`; if symbols are missing at the
  final link of `gchp`, make it `PUBLIC`.
- Tapenade C compile uses `-fPIE` (GNU/Intel ok).

Keep compile fixes minimal and list each in `ADJOINT_PORT_NOTES.md`.

## 4. Run directory

Create a 14.8.0 carbon run directory with the 14.8.0 `createRunDir.sh` (do not reuse 14.5.3
config files: 14.8.0 changed `GCHP.rc`, `ExtData.rc`, `HISTORY.rc`, `geoschem_config.yml`).
Then add the adjoint keys from the 14.5.3 run directories (same names and defaults):
`MODEL_PHASE`, `REVERSE_TIME`, `ADV_ADJ_METHOD: reverse_wind_pure`, `ADV_ADJ_WEIGHT: dpadv`,
`ADV_ADJ_MASSFLUX: forward`, `ADV_ADJ_RESCALE: 1`, `FD_TYPE`, `FD_SPEC`, `FD_STEP`, `IFD/JFD`
or `FD_LAT/FD_LON`, `LFD`, `FD_LAT/LON_MIN/MAX`, `ADJOINT_VAL`, `ADJOINT_FROM_FILE`,
`OCO2_FORCING_DIR`, `ADJ_HEMCO_SFLUX_*`, and the `SpeciesAdj`/`SurfaceFluxAdj` History
collections. Use the forward `ExtData.rc` unchanged for the adjoint (MAPL fix ported).

14.8.0-specific settings to check:
- 14.8.0 replaced `METEOROLOGY_VERTICAL_INDEX_IS_TOP_DOWN` with `USE_EXTDATA2G`,
  `MET_WIND_IS_TOP_DOWN`, `MET_HUMIDITY_IS_TOP_DOWN`, `MET_NONADVECTION_IS_TOP_DOWN`; keep the
  values `createRunDir.sh` writes.
- **Restart `DELP_DRY`:** 14.8.0 scales all tracers on the first AdvCore call by
  restart dp / met dp if the restart has a non-zero `DELP_DRY`. The transpose of this is not
  implemented (notes 2.5). For all FD / adjoint tests use restarts **without** `DELP_DRY`
  (or with it zeroed) unless you are testing this point on purpose; record which you used.

## 5. Rules

- Push only to the `kaysuselj/*` forks, branch `adjoint_port_14.8.0`; never to
  `geoschem/*` or `GEOS-ESM/*`; no force-pushes; ask the user before every push.
- When you change a submodule, commit there first, then update the pointer in GCHP.
- Do not "clean up" adjoint numerics (method, weighting, rescale, debug blocks) in this pass.
- The forward model must stay identical to standard GCHP 14.8.0 (apart from the first-step
  advection fix). Never change forward behaviour to make the adjoint easier.
- Judge the port by adjoint vs finite differences in 14.8.0, against the 14.5.3 numbers; do not
  expect bit-identity with 14.5.3 (upstream science changed, e.g. post-advection T/SPHU).
- Pitfalls from 14.5.3 still apply: vertical order (restarts/checkpoints top first, History
  surface first, `LFD` surface based), `SpeciesAdj_CO2` per kg/kg (x1.5195 vs mol/mol),
  end-time adjoint record is zero, checkpoints are consumed by adjoint runs.

## 6. Tests (in this order; later ones need the earlier ones)

| # | Test | Pass criterion |
|---|---|---|
| T0 | Forward run, ADJOINT build vs plain build, a few days, advection on | Bit-identical from step 1 (checks the first-step fix in GCHPctmEnv and AdvCore, and that nothing adjoint-only leaked into the forward). Most important regression test of the port. |
| T1 | Forward vs adjoint met (the 17 fields of 14.5.3 plus `PS2`, `SPHU2`, `TMPU2`, `DELP_DRY`) and emissions, step by step | Bit-identical with the forward `ExtData.rc` |
| T2 | Forcing loader on 6 vs 24 cores | Bit-identical; every observation deposited once |
| P1a | Column tests, advection off: initial-CO2 adjoint vs FD | 14.5.3: max deviation 1.2e-6. Also exercises the AIRQNT adjoint guard (notes 2.1) and DELP_DRY restore (2.2) |
| P1b | Surface-flux adjoint vs FD, PBL + convection | 14.5.3: max deviation 1.2e-5, all columns (convection, notes 2.3) |
| T6 | One-step advection matrix, all cells | Row sums 1 to 1e-5 (median), diagonal error ~1e-5 |
| P2 | Idealised winds, 5 days, inside a face | 0.7-2% |
| P2c | Polar-face crossing | 10.7% |
| P2a | Seam / corner, flat ocean, 5 days | 6% seam, 19-32% corner |
| P2e | Real winds, 1 month, block receptors | 3-13% (point receptors 14-36%) |
| P3 | One-month OCO-2 cost, full physics | Slope 0.98 (flux), 0.99 (initial CO2) |

P2e and P3 reference numbers are from the paper binary; the user is rerunning the advection
tests with `gchp.ref_20261014` and may pass on updated numbers. If a 14.8.0 result is clearly
worse than 14.5.3, find the first test in the table that degrades and look at the matching item
in `ADJOINT_PORT_NOTES.md` section 2 before changing anything.

## 7. Open items you may be asked to do after the tests

- Transpose of the 14.8.0 first-step `DELPDRY` scaling (notes 2.5).
- OSSE sampling with the MAPL 2.69 trajectory sampler instead of `track_file` (notes 3).
- Fix the `NA = NFD` species/advected-index mix in the convection adjoint if species other
  than CO2-as-species-1 are needed (notes 2.3).
