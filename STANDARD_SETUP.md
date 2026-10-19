# GCHP adjoint / CO2 4D-Var: standard setup (14.8.0, optimized)

Decision 2026-10-16: **all new runs use GCHP 14.8.0 with the optimized (reduced-checkpoint) adjoint.**
Decision 2026-10-19: **14.5.3 is stale.** It is kept unchanged in `/nobackupp27/ksuselj1/gchp_14.5.3_adjoint_surfaceF/`
only so we can go back and check if 14.8.0 gives strange results. Everything 14.8.0 lives in
`/nobackupp27/ksuselj1/gchp_14.8.0_adj/` and does not read from the 14.5.3 tree.

## 0. Where things are (gchp_14.8.0_adj/)

| path | contents |
|---|---|
| `code/` | 14.8.0 adjoint source (git, kaysuselj/GCHP `adjoint_port_14.8.0`), builds `build_adj/`, `build_plain/` |
| `gchp_4Dvar/` | 4D-Var driver and Python tools (git clone of kaysuselj/gchp_4Dvar) |
| `4dvar_osse/setup_exp.sh` | **the script to create an OSSE experiment** |
| `4dvar_osse/exps/<name>/` | **where each OSSE runs**: `4dvar_optimizer.osse.<name>.run` (submit with qsub), `forward_run/`, `adjoint/`, `logs/`, `plots_osse/`, `4dvar_output_osse/` (state, resume scripts) |
| `4dvar_osse/control_files_148*/`, `truth148/`, `obs148/`, `fd_check/`, `make_ctrl148.sh` | control files, truth run, observations, gradient check |
| `inputs/` | OSSE surface fluxes, idealized met, C24 base restart + OCO-2 track (moved/copied out of the 14.5.3 tree) |
| `tests_final/`, `tests_final_c48/` | adjoint tests at C24 and C48 (14.5.3 reference inputs copied to `tests_final/common/ref153/`) |
| `summary/` | slides and figures |

The old 14.5.3 paths of moved items are symlinks to the new places, so the stale 14.5.3 scripts still run.

## 1. The standard

| item | standard |
|---|---|
| Model | GCHP 14.8.0 with the adjoint port (MAPL 2.69, ESMF 8.6.1 from source) |
| Adjoint binary | `code/build_adj/bin/gchp.ckpt_20261016` |
| Forward-only binary (FD checks, truth runs that need no adjoint checkpoints) | `code/build_plain/bin/gchp.plain_20261014` (same forward) |
| Checkpointing | `ADJ_CKPT_REDUCED: 1`. The adjoint reads only the fields it needs. The gradient is bit-identical to full checkpoints, 40–45 % faster. |
| Advection adjoint | `ADV_ADJ_METHOD: reverse_wind_pure`, `ADV_ADJ_WEIGHT: airmass`, `ADV_ADJ_MASSFLUX: forward` |
| Resolution / step | C24, 1800 s (C48 being tested, see §4) |
| 4D-Var nodes | 2 Broadwell nodes (48 cores) |
| 4D-Var stopping | `--ftol 0.02` (stop when J improves by < 2 %; E1.3 gets there around iteration 10) |
| Plots | `--plot-every 5` (plus iteration 1, the convergence plot and the end) |
| Restart | `DELP_DRY = 0` (14.8.0 first-step pressure scaling has no adjoint) |

Expected cost of a 3-month C24 OSSE: about 1.9 h per iteration and about 11 iterations, so **about 21 h**.
The 14.5.3 full-checkpoint setup took about 53 h.

## 2. How to start new runs

- **4D-Var OSSE:**
  `4dvar_osse/setup_exp.sh <name> --prior-err "<opts>"`
  - The defaults are now the standard: `--gchp 14.8`, which uses `control_files_148_ckpt_reduced/` and the
    14.8.0 truth observations `obs148/`, plus `--ftol 0.02`, `--plot-every 5` and `--nodes 2`. The old
    script is kept as `setup_exp.sh.pre_default148_20261016`.
  - `--gchp 14.8-fullckpt` gives full checkpoints (only for comparisons).
  - `--res 48` (2026-10-19) sets up a C48 experiment: `control_files_148_c48/` (MERRA-2 2x2.5 met, reduced checkpoints),
    C48 truth observations `obs148_c48/`, E on the C48 grid, longer driver time budgets. Default `--res 24`.
  - 14.5.3 is no longer available from this script (stale tree: `gchp_14.5.3_adjoint_surfaceF/4dvar_final/setup_exp.sh`).
  - Then: `cd exps/<name>; qsub 4dvar_optimizer.osse.<name>.run`.
- **C24 adjoint tests** (`tests_final/`): the drivers (`run_P1.pbs`, `run_P2.pbs`, `run_ADV.pbs`,
  `run_P2lin.pbs`) now default to `gchp.ckpt_20261016`. The tests read full checkpoints, so they leave
  `ADJ_CKPT_REDUCED` off; the forward and adjoint are the same as before.
- **C48 tests** (`tests_final_c48/`): already on the optimized binary.

## 3. Why 14.8.0 (validation, all against finite differences)

| test | 14.8.0 result |
|---|---|
| P1a/b columns (PBL, + convection) | exact (~1e-6); with the OSSE fluxes the same response as 14.5.3 |
| P1c/d (decomposition, timing) | bit-identical |
| P2a, P2c, P2d, P2f advection | identical to 14.5.3 |
| T2 / P2e blocks, real winds | 4.4–5.5 % / 4–13 % (same as 14.5.3) |
| OSSE FD gradient check, 12 cells | slope 0.994, weighted error 1.07 % |
| OSSE iteration-1 gradient vs 14.5.3 | correlation 1.000000 |
| E1.3 4D-Var J vs 14.5.3 | ≤ 0.13 % different through iteration 4 |
| Reduced checkpoints | gradient bit-identical, 40–45 % faster |
| C48 columns (P1) | exact (1.4e-6 / 1.2e-5) |

Details: `code/ADJOINT_PORT_NOTES.md` (§5), `4dvar_osse/SPEEDUP_REPORT.md` (§6).

## 4. Runs not on the standard (left as they are)

- `gchp_14.5.3_adjoint_surfaceF/4dvar_final/exps/E1.3_const_sat` (14.5.3, finished) and `4dvar_osse/exps/E1.3_const_sat_148` (14.8.0 with full
  checkpoints). Both were started before the switch and keep their own copies. They are the
  14.5.3-vs-14.8.0 comparison.
  - `E1.3_const_sat_148` resumes with `RESTART=true` when its walltime runs out. It may be resumed with
    the optimized binary, since the forward is identical.

## 5. Known limits and next steps

- The advection adjoint is a reversed-wind approximation. Errors are concentrated at cube-face edges and
  corners (blocks 4–13 %). An exact (discrete) adjoint is scoped in `FULL_ADVECTION_ADJOINT_SCOPE.md`,
  option B first, from the week of 2026-10-19. The devkit is `adv_adjoint_devkit_C24.tar.gz`.
- C48 and then C90: `tests_final_c48/PLAN_C48_C90.md` (status marks, updated per step).
- Later on 14.8.0: spatially correlated background error (B); bilinear observation sampling as an option
  next to nearest-neighbour; per-rank in-memory checkpoints and Aitken nodes (needed for C90).
- Nothing from the 14.8.0 port or the optimization is committed yet. Commits and pushes go to the
  kaysuselj forks only, after approval.

All paths are relative to `/nobackupp27/ksuselj1/gchp_14.8.0_adj/` (= `/home5/ksuselj1/nobackup/gchp_14.8.0_adj/`).
