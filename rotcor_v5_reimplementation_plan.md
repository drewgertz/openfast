# RotCor v5 Reimplementation Plan

## Goal
Reimplement RotCor natively for OpenFAST v5 with minimal coupling risk, then add complexity in validated steps.

## Scope
- AeroDyn only.
- BEM-based path only for active RotCor.
- Preserve default behavior when `RotCor = 0`.

## Branch
- Working branch: `feat/rotcor-v5-reimpl`

## Phase 0: Baseline and Guardrails
- Build OpenFAST on this branch with no new logic changes.
- Confirm a baseline case with `RotCor = 0` reproduces upstream-v5 behavior.
- Record baseline channels for comparison:
  - Blade node aerodynamic loads/coefs (Cl/Cd-related channels)
  - Rotor power/thrust channels

## Phase 1: Input and Runtime Wiring (No Physics Change)
- Files:
  - `modules/aerodyn/src/AeroDyn_IO.f90`
  - `modules/aerodyn/src/AeroDyn.f90`
  - `modules/aerodyn/src/AirfoilInfo_Types.f90`
- Tasks:
  - Keep `RotCor` parse/validation clear and explicit.
  - Ensure runtime state for RotCor is initialized and updated only on BEM path.
  - Add one-time diagnostic message in initialization to report effective RotCor mode.

## Phase 2: Minimal Physics Path (RotCor=1 Only)
- Files:
  - `modules/aerodyn/src/BEMTUncoupled.f90`
  - `modules/aerodyn/src/AirfoilInfo.f90`
  - `modules/aerodyn/src/AirfoilInfo_RotCor.f90`
- Tasks:
  - Implement and validate only Snel (`RotCor=1`) first.
  - Keep this as a direct runtime path, not an aggressive precompute path.
  - Ensure local snel factor is updated each call on BEM path.

## Phase 3: Expand Models
- Add and validate:
  - `RotCor=2` (Schepers)
  - `RotCor=3` (Snel + Gertz Cd)

## Phase 4: Optional Precompute Optimization
- Only after physics validation:
  - Add precompute table optimization with strict parity tests versus runtime path.

## Validation Matrix
- Cases:
  - `Wake_Mod=1`, `RotCor=0`
  - `Wake_Mod=1`, `RotCor=1`
  - `Wake_Mod=1`, `RotCor=2`
  - `Wake_Mod=1`, `RotCor=3`
  - `Wake_Mod=0` and `Wake_Mod=3` with `RotCor>0` (must stay effectively inactive)
- Checks:
  - `RotCor=0` remains unchanged vs baseline.
  - `RotCor>0` produces non-identical aerodynamic channels.
  - No segfaults on constant/non-constant airfoil tables.

## Exit Criteria
- OpenFAST builds cleanly.
- `RotCor=0` behavior preserved.
- `RotCor>0` has measurable, physically plausible impact in BEM cases.
- No regression crashes in prior failing scenarios.

## Notes
- Do not optimize early.
- Keep each phase mergeable and testable independently.
