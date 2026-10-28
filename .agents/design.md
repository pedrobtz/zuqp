# RFC 0011: `zuqp`, dense quadratic programming in R with DAQP

- **Status:** Adopted. This file is the package's `.agents/design.md`;
  `.agents/roadmap.md` carries the stages and their exit criteria. Each
  decision in section 16 stays the maintainer's until a stage records it.
  Revised 2026-10-13 against the DAQP source (see section 21).
- **Date:** 2026-10-09; revised 2026-10-13.
- **Package:** `zuqp` (repository skeleton created; neither `zuqp` nor
  `daqp` was in CRAN's package index on 2026-10-13). Vendors DAQP.
- **Source:** the DAQP repository (`darnstrom/daqp`) at `master` commit
  `361ce54` (2026-10-08; "Eliminate the equality constraints of a problem
  with a factored Hessian"); `master` had moved to `ce6204e` (2026-10-10,
  one commit later) when this revision was made. The latest tag is
  `v0.10.3` (`29a8987`, 2026-10-02), five commits behind `361ce54`; those
  five commits touch seven files in `src/` and `include/` and add the
  exit flag `DAQP_EXIT_FEASIBLE` (§5), so the tag and the commit are not
  interchangeable (§9). Read: `include/api.h`, `include/types.h`,
  `include/constants.h`, `include/utils.h`, `src/*.c` (ten files; 6,438
  lines with the thirteen headers), `docs/index.md`,
  `docs/docs/{settings,c,miqp}.md`, `LICENSE`, `CMakeLists.txt`, and the
  Python (`interfaces/daqp-python/daqp.pyx`, `setup.py`) and Julia
  (`interfaces/daqp-julia/src/api.jl`) interfaces; Arnström, Bemporad
  and Axehill, "A Dual Active-Set Solver for Embedded Quadratic
  Programming Using Recursive LDLᵀ Updates", IEEE TAC 67(8), 2022; the
  CRAN Optimization task view, which lists no DAQP binding; the
  alternatives review that proposed this package ("DAQP for quadprog: a
  small, library-free C solver; dense active-set, warm starts"). The
  reading was of source, not of a running build; section 13 says what
  the first stage measures.

## 1. What zuqp is

`zuqp` solves convex quadratic programs

```text
minimise    ½ xᵀ H x + fᵀ x
subject to  l  ≤ x  ≤ u
            bl ≤ A x ≤ bu
```

through DAQP, a dual active-set solver in plain C with no dependencies,
written for the small dense problems that model predictive control
solves thousands of times a second. In R those are the problems
`quadprog::solve.QP()` has served since the 1990s: portfolio
optimisation, index tracking, constrained least squares, hedging, small
risk-budgeting and calibration steps run in loops.

It is a *solver* package, like `zusat` and `zusmt`. Beyond plain QPs it
exposes what DAQP adds over `quadprog`: box bounds handled as bounds,
two-sided constraints, equality rows, a positive semidefinite `H` (so
linear programs too) through proximal-point iterations, warm starts for
sequences of nearby problems, and binary constraints for small
mixed-integer QPs.

The one-line statement that governs every decision below:

> **Solve small dense QPs exactly as quadprog users expect, faster, and
> make the repeated solve, not the single solve, the fast path.**

## 2. Scope

**In v0.1.0:**

- `qp_solve()`: one QP, with bounds, two-sided general constraints and
  equalities (`bl == bu`).
- Positive semidefinite `H`, including `H = NULL` (linear programs),
  through DAQP's proximal-point iterations.
- `qp_model()`: a reusable model whose `f`, bounds and constraint
  vectors (and optionally `H` and `A`) are updated in place and
  re-solved warm.
- Binary constraints (`binary =`), solved by DAQP's branch and bound.
  DAQP's branch and bound runs only on the plain active-set path: with a
  singular `H` or an LP it returns `DAQP_EXIT_NONCONVEX`
  (`src/api.c`, `daqp_solve`), so binary constraints need a positive
  definite `H` (§6).
- A `quadprog`-compatible front end, `qp_solve_quadprog()`, taking
  `solve.QP()`'s `Dmat, dvec, Amat, bvec, meq` and returning its list
  shape, so migration is a rename.
- Exit status as a classed result field, never a silent wrong answer.
- Every field of `DAQPSettings` reachable through `qp_settings()`,
  including the branch-and-bound pruning tolerances and the
  equality-reduction policy (§5).

**Later, when asked:**

- Hierarchical QPs and affine variational inequalities (DAQP has both:
  `setup_daqp_hiqp()`, `daqp_avi()`).
- Soft constraints with weights (`daqp_set_soft_weights()`).
- Primal and dual starting points for one-shot solves
  (`daqp_primal_init_active()`, `daqp_dual_init_active()`).
- An `ROI` plugin (`ROI.plugin.zuqp`), so ROI users reach DAQP.
- Batch solving of many independent problems, possibly in parallel.

**Never:**

- Sparse problems. DAQP is dense by design and its documentation sends
  large sparse problems to OSQP (`docs/index.md`); so does zuqp
  (`osqp`, `piqp` and `clarabel` are on CRAN).
- A modelling language. `CVXR` and `ROI` exist.

## 3. Position in the `zu*` family

The rows for `zu-family.md` (alignment R1):

| | zuqp |
| --- | --- |
| Role | solver |
| R prefix | `qp_` |
| Info function | `zuqp_info()` |
| Root condition class | `zuqp_error` |
| Public C prefix | none; internal `zqp_` / `ZQP_` |
| From C | no C API in 0.1.0 |
| Hides symbols | yes, from Stage 0 (R3) |
| Vendored code | DAQP (MIT) |
| `Depends: R` | 4.1 (R2) |
| Language | en-GB (R2) |

The solver packages are the family's weakest on conventions (alignment
R9: `zusat` and `zusmt` raise bare `stop()` and have no info function).
zuqp starts with R9 satisfied, so it can serve as their template.

**Relationships.** None with other family packages: DAQP needs no
buffers, parsing or big integers. `zuprima` (RFC 0012) is the sibling for
derivative-free problems, and both follow the same result shape (§5).

**Consumers.** Portfolio and risk code (mean-variance, minimum tracking
error, risk parity steps), constrained regression, MPC-style
simulations in R, and every package that calls `solve.QP()` in a loop
(about two hundred reverse dependencies of `quadprog`).

**What exists elsewhere.**

- `quadprog` (CRAN; last release 2019): Goldfarb–Idnani in Fortran
  from the 1980s. Strictly positive definite `D` only; no bounds (they
  must be written as constraint rows); one-sided inequalities; no warm
  start.
- `osqp`, `piqp`, `clarabel` (CRAN): ADMM and interior-point methods,
  strong on large or sparse problems, but not warm-startable in the
  active-set sense and slower per solve on small dense ones.
- `ROI.plugin.qpoases` (CRAN): qpOASES, an online active-set solver in
  C++ reached through the ROI framework. The closest in kind; DAQP
  differs in being plain C with no dependencies, in its binary
  constraints, and in the speed its authors report on small problems.
- `kernlab::ipop()`, `LowRankQP`, `limSolve`: older interior-point and
  special-purpose routines.

The gap: a small, fast, dependency-free dense QP solver with warm starts
and a drop-in path for `quadprog` users.

## 4. Architecture

```text
R (H, f, A, bounds, settings)
   |
   v
validate and pack  (R/qp-solve.R, src/zqp_pack.c)
   - shapes, finiteness, symmetry of H to a tolerance
   - bounds and general constraints concatenated as DAQP expects:
     blower = [l; bl], bupper = [u; bu], ms = n always (§6)
   - A transposed into row-major order (DAQP's layout; the Python
     interface passes a C-contiguous array and the Julia interface
     assigns A' for the same reason)
   - sense flags per constraint: ACTIVE + IMMUTABLE for an equality,
     BINARY for a binary row, 0 otherwise
   |
   v
DAQP  (src/daqp/, vendored, compiled with -DPROFILING)
   - setup_daqp_main() + daqp_solve() for one solve (what
     daqp_quadprog() does, inlined so the workspace is readable)
   - setup_daqp() + daqp_update_ldp() + daqp_solve() for models
   |
   v
result  (src/zqp_result.c)
   - x, objective, duals (lam), iterations, nodes, times, whether
     proximal regularisation was used (work->n_prox), exit flag
     mapped to a status string
```

A model (`qp_model()`) is an external pointer to a `DAQPWorkspace` plus
the problem arrays it points into, all allocated in C and freed by one
finaliser. R vectors are copied in, never aliased, so R's copy-on-modify
cannot change a workspace behind DAQP's back.

One-shot solves do not call `daqp_quadprog()`. That function builds a
workspace on its stack, solves, and frees it, so nothing but
`DAQPResult` survives; zuqp does the same three steps itself
(`setup_daqp_main()` with `DAQP_UPDATE_unconstrained |
DAQP_UPDATE_eliminate`, `daqp_solve()`, the two `free_daqp_*()` calls)
so that it can read `work->n_prox` before freeing and report whether
the problem was regularised (§8). `daqp_update_ldp()` is declared in
`include/utils.h`, not `include/api.h`; the Python and Julia interfaces
call it the same way, so it is API in practice (§9).

DAQP calls no R API and prints nothing: `src/api.c` and `src/utils.c`
include `<stdio.h>` but call nothing from it, and there is no `printf`,
`exit`, `abort` or `rand` in `src/` (checked at `361ce54`). The timer
exists only under `#ifdef PROFILING` (`src/utils.c`, `src/api.c`): it is
`clock_gettime(CLOCK_MONOTONIC)` on POSIX and `QueryPerformanceCounter`
under `_WIN32`. Without `PROFILING`, `time_limit` is never checked and
the result's `solve_time` and `setup_time` are zero. DAQP's own CMake
build and the Python package both define it; so does zuqp (D10, §13).

## 5. Public R API (complete v0.1.0 surface)

```r
# one problem
qp_solve(H, f, A = NULL, bl = NULL, bu = NULL,
         lower = NULL, upper = NULL, binary = NULL,
         settings = qp_settings())

# quadprog's interface, quadprog's result shape
qp_solve_quadprog(Dmat, dvec, Amat, bvec, meq = 0)

# repeated problems
m <- qp_model(H, f, A, bl, bu, lower, upper, binary = NULL,
              settings = qp_settings())
qp_update(m, f = , bl = , bu = , lower = , upper = , H = , A = )
qp_solve_model(m)

qp_settings(primal_tol = 1e-6, dual_tol = 1e-12, iter_limit = 10000,
            time_limit = 0, eps_prox = -1e-6, eq_reduction = 0,
            rel_subopt = 0, abs_subopt = 0, fval_bound = 1e30, ...)
zuqp_info()
```

`qp_settings()` takes every field of `DAQPSettings` under DAQP's name
with DAQP's default (`daqp_default_settings()`), so DAQP's settings page
documents zuqp's; the signature above shows the ones users will reach
for. Two defaults differ from what a `quadprog` user might expect and
are documented up front: `eps_prox = -1e-6` means *automatic* proximal
regularisation of a singular `H` (a positive value forces it, `0`
disables it), and `time_limit = 0` means *no* limit (`Inf` is accepted
and mapped to `0`).

**Result** (a list of class `zuqp_result`, the shape `zuprima` shares):

| Field | Meaning |
| --- | --- |
| `x` | solution |
| `value` | objective at `x` |
| `status` | one of the statuses in the table below |
| `success` | `status %in% c("optimal", "optimal_inexact")` |
| `duals` | multipliers, bounds first then rows, signed by side (Stage 1 records the sign convention from `ldp2qp_solution`) |
| `active` | logical, which constraints have a nonzero multiplier (weakly active constraints read `FALSE`; documented) |
| `regularised` | `TRUE` if DAQP used proximal-point iterations (`work->n_prox > 0`) |
| `iterations`, `nodes` | DAQP's counters |
| `solve_time`, `setup_time` | seconds, from DAQP's timer |
| `code` | DAQP's raw exit flag |

Statuses, one per exit flag in `include/constants.h` at `361ce54`:

| `code` | DAQP flag | `status` |
| --- | --- | --- |
| 5 | `FEASIBLE` | `feasible` (branch and bound stopped by a limit with an integer-feasible incumbent; not proven optimal) |
| 4 | `OPTIMAL_INEXACT` | `optimal_inexact` (a constraint is violated by more than `primal_tol` after cycling) |
| 2 | `SOFT_OPTIMAL` | `soft_optimal` (reserved; needs soft constraints) |
| 1 | `OPTIMAL` | `optimal` |
| −1 | `INFEASIBLE` | `infeasible` |
| −2 | `CYCLE` | `cycling` |
| −3 | `UNBOUNDED` | `unbounded` |
| −4 | `ITERLIMIT` | `iteration_limit` |
| −5 | `NONCONVEX` | `nonconvex` (also binary constraints with a singular `H`) |
| −6 | `OVERDETERMINED_INITIAL` | from setup: a `zuqp_setup_error`; from a solve: `overdetermined_initial` |
| −7 | `TIMELIMIT` | `time_limit` |
| −8 | `UNSUPPORTED` | `unsupported` (hierarchies with prox or AVI; cannot occur in 0.1.0) |

An unknown flag (a re-vendor that adds one) maps to `unknown` with the
raw `code`, never to an error, so the package degrades rather than
breaks.

`qp_solve()` never errors for an infeasible or unbounded problem: the
status says so. It errors (classed) only for invalid input and for a
setup DAQP refuses (§10).

`qp_solve_quadprog()` maps `solve.QP()` exactly: `Amat` columns are
constraints, the first `meq` are equalities, `bvec` is the lower side,
and the result has `solution`, `value`, `unconstrained.solution`,
`iterations`, `Lagrangian` and `iact`, as `quadprog` returns them. Two
fields have no DAQP counterpart and are documented: `iterations` is
`quadprog`'s length-two vector (`c(total, dropped)`), returned as
`c(iterations, NA_integer_)`; `iact` lists the active constraints in
ascending index order, not in the order `quadprog` added them.
`Lagrangian` is non-negative, converted from DAQP's signed multipliers.
On an infeasible problem it raises `zuqp_infeasible`, as `quadprog`
stops.

## 6. Problem mapping

| User input | DAQP |
| --- | --- |
| `lower`, `upper` (`±Inf` allowed, `NULL` means unbounded) | simple bounds: `ms = n`, rows `1…n` of `blower`/`bupper`, `±DAQP_INF` where unbounded |
| `A`, `bl`, `bu` | general rows `n+1…m`; `±Inf` sides become `±DAQP_INF` |
| `bl[i] == bu[i]` | an equality: `sense = DAQP_ACTIVE + DAQP_IMMUTABLE` (5), set explicitly as DAQP's C docs do |
| `binary = idx` (bound or row indices) | `sense = DAQP_BINARY` (16): `bl` or `bu`, nothing between; at most `n` of them (`setup_daqp_bnb`) |
| `H` positive semidefinite | DAQP's automatic regularisation (`eps_prox < 0`), reported as `regularised = TRUE` |
| `H = NULL` | linear program: `setup_daqp_ldp()` marks all `n` directions proximal; `f` must then be given |

`ms` is always `n`: DAQP's first `ms` constraints are bounds on the first
`ms` *variables* (`blower[1:ms] ≤ x[1:ms] ≤ bupper[1:ms]`, Julia docs),
so a partial bound set cannot be packed as a shorter block without
reordering variables. Unbounded variables get `±DAQP_INF` instead
(D11). A variable bound on one side only is therefore free, not an
error.

Infinity: DAQP's constant is `DAQP_INF = 1e30`; `±Inf` from R is mapped
to it, and a finite bound at or beyond it is a `zuqp_invalid_argument`
rather than a silent clamp.

Binary constraints are DAQP's general form `A x ∈ {bl, bu}`; `x_i ∈
{0, 1}` is the bound case with `lower = 0, upper = 1` and is the
documented example (open question 3).

## 7. Warm starts and models

`qp_model()` exists because repeated solves are where an active-set
method wins: after `qp_update()` changes `f` or the bounds, DAQP restarts
from the previous active set and usually needs a handful of iterations.

`qp_update()` calls `daqp_update_ldp()` with the mask of what changed,
as the Python and Julia interfaces do: `DAQP_UPDATE_v` for `f`,
`DAQP_UPDATE_d` for any bound vector, `DAQP_UPDATE_Rinv` for `H`,
`DAQP_UPDATE_M` for `A`. Changing `H` or `A` refactorises the working
set and is documented as the expensive case. In that case DAQP also
requires `DAQP_UPDATE_sense`, and keeps the previous active set only if
`qp->sense` *is* `work->sense` (`docs/docs/c.md`: "point the problem at
the workspace sense before updating"); a separate sense array replaces
the state. The model therefore sets `qp.sense = work.sense` once after
setup and never passes a separate sense array (D12).

Consequences, all documented:

- A row's equality status is fixed for the life of a model: an update
  that turns `bl[i] == bu[i]` into `bl[i] < bu[i]`, or the reverse, is
  a `zuqp_invalid_argument` ("rebuild the model"). Keeping the sense
  array stable is what makes the warm start valid.
- A model is never equality-reduced: DAQP's automatic reduction
  (`eq_reduction = 0`) applies only to solves that start from scratch,
  because warm solves are too short to recover its cost
  (`docs/docs/c.md`). `qp_solve()` gets the reduction through
  `DAQP_UPDATE_eliminate`; `qp_settings(eq_reduction = 1)` forces it on
  a model for callers who update rarely and solve from scratch.
- `H` and `A` must keep their dimensions; `qp_update()` checks before
  touching the workspace, so a failed update leaves the model solvable.

Typical loop:

```r
m <- qp_model(Sigma, -mu, A = matrix(1, 1, n), bl = 1, bu = 1,
              lower = rep(0, n), upper = rep(0.1, n))
for (t in seq_len(T)) {
  qp_update(m, f = -mu_t[t, ])
  w[t, ] <- qp_solve_model(m)$x
}
```

## 8. Numerical policy

- `H` is symmetrised (`(H + t(H)) / 2`) after a check that it was
  symmetric to `1e-10` relative; otherwise a `zuqp_invalid_argument`.
- Non-finite entries in `H`, `f` or `A` are errors; `±Inf` only in bounds.
- Positive semidefinite problems: DAQP's default `eps_prox = -1e-6`
  regularises a singular factorisation automatically, so zuqp has no
  retry of its own. It reports what happened (`regularised = TRUE`) and
  never silently changes the problem. With `eps_prox = 0` a singular
  `H` is DAQP's `nonconvex` status (Stage 4 confirms the flag).
- Tolerances are DAQP's defaults, exposed in `qp_settings()` with
  DAQP's names (§5).

## 9. The vendored library, as zuqp uses it

- `src/daqp/`: the ten `src/*.c` files and thirteen `include/*.h` at a
  pinned commit, plus `LICENSE` into `inst/COPYRIGHTS`. Neither
  directory's `CMakeLists.txt` is copied.
- The pin is a commit, not a tag (D13): `v0.10.3` predates
  `DAQP_EXIT_FEASIBLE` and the equality-elimination fixes of the five
  commits that follow it. Stage 0 pins the first tag that includes
  `361ce54` if one exists by then, otherwise `361ce54` or a later
  `master` commit that Stage 0 has built; `src/daqp/VERSION` records
  the commit (and the tag when there is one) and the SHA-256 per file.
- Built with `DAQP_SINGLE_PRECISION` undefined (double) and
  `-DPROFILING` defined (D10); `DAQP_SOFT_WEIGHTS` is left at DAQP's
  default (it only adds workspace fields that are allocated on demand).
  The `codegen/` and `interfaces/` directories are not vendored.
- `tools/vendor-daqp` fetches a commit, copies the files, and writes
  `src/daqp/VERSION`.
- Patches: none expected. DAQP has no stdio calls and no `exit`; Stage 0
  confirms with CRAN's compiled-code check and records any warning.
- All DAQP calls go through `src/zqp_daqp.h`, which includes `api.h` and
  `utils.h` (for `daqp_update_ldp()`), so a re-vendor that renames a
  function touches one file.

## 10. Errors

Every condition inherits `zuqp_error`: `zuqp_invalid_argument`
(shapes, non-finite data, asymmetric `H`, `bl > bu`, a finite bound
beyond `1e30`, more than `n` binary constraints, `H = NULL` without `f`,
an equality row whose status an update would change, binary
constraints with `H = NULL`), `zuqp_infeasible` (only from
`qp_solve_quadprog()`, to match `quadprog`), `zuqp_setup_error` (DAQP's
setup returned a negative flag: `OVERDETERMINED_INITIAL` from
`setup_daqp_bnb()` or `daqp_update_ldp()`, `NONCONVEX` from the
factorisation with `eps_prox = 0`; the message carries the status
name), and `zuqp_internal_error`. Solve outcomes are statuses, not
conditions (§5).

## 11. Limits

- Dimensions: DAQP uses `int`; `n` and `m` above `2^31 − 1` are refused.
  In practice the solver is meant for `n` in the tens to low thousands;
  `qp_solve()` warns (`zuqp_large_problem`) above `n = 5000`, pointing to
  `osqp`.
- Iterations: `iter_limit` defaults to DAQP's 10,000.
- Time: `time_limit` passes through and depends on `PROFILING` (§4);
  the CI matrix tests that it terminates a solve on every platform.
  Interrupts are checked between solves in batch code, never inside
  DAQP, which has no interrupt hook; a long branch and bound is stopped
  with `time_limit`, not `Ctrl-C`.
- Branch and bound: bounded by `iter_limit` and `time_limit`; DAQP has no
  separate node limit at `361ce54`. `rel_subopt`, `abs_subopt` and
  `fval_bound` prune the tree and are the documented way to trade
  optimality for time; a limit reached with an incumbent returns
  `feasible`, not an error.
- Memory: DAQP allocates with `malloc` and does not check for failure.
  zuqp validates sizes before setup; an allocation failure inside DAQP
  is out of scope, as it is for every vendored C library in the family.

## 12. Memory model

One external pointer per model owns the workspace and every array DAQP
reads; the finaliser calls `free_daqp_ldp()` and then
`free_daqp_workspace()`, which itself frees the branch-and-bound tree
(`free_daqp_bnb()`) and, when `work->settings` is not `NULL`, the
settings struct. zuqp owns its settings struct, so it sets
`work->settings = NULL` before that call, exactly as `daqp_quadprog()`
does; forgetting this is a double free. One-shot `qp_solve()` builds
the same structure on the C stack and heap and frees it before
returning, with R allocation of the result done after DAQP returns, so
an R error cannot leave DAQP memory behind.

## 13. Build, portability and CRAN

- C99 under the family's lint flags; DAQP's sources compile under the
  same flags or get a recorded, narrow exception.
- `DESCRIPTION`: no `Imports`, no `LinkingTo`; `Suggests: quadprog,
  testthat, withr` (`quadprog` as the oracle).
- `src/Makevars`: `PKG_CFLAGS = $(C_VISIBILITY)`, `PKG_CPPFLAGS =
  -DPROFILING -I daqp/include`; `R_init_zuqp` the only export; no BLAS
  or LAPACK needed. On Windows `include/utils.h` pulls in `<windows.h>`
  for the timer; Stage 0 confirms it builds under Rtools and that
  `clock_gettime` needs no extra library on the Linux and macOS legs.
- Stage 0 measures warnings under CRAN's gcc and clang flags and the
  compiled-code NOTE before anything else.

## 14. Testing

- **The oracle is `quadprog`**: for every random strictly convex problem
  with a known feasible point, `qp_solve()` and `solve.QP()` agree on `x`
  and the objective to `1e-8` relative, and `qp_solve_quadprog()`
  returns the same list.
- **Maros–Mészáros** small dense problems (the standard QP test set)
  against published optimal objectives, for the subset DAQP's
  documentation targets.
- **Warm starts**: a model updated `k` times gives the same `x` as `k`
  cold solves, with fewer iterations (recorded, not asserted); an
  update of `H` keeps the active set (the sense-aliasing contract of
  §7) and an update that changes a row's equality status is refused.
- **Edge cases**, each a permanent regression: infeasible bounds,
  unbounded LP, `H = NULL` with bounds, singular `H` (`regularised`
  reads `TRUE`), duplicate constraint rows, all-equality problems,
  `n = 1`, `±Inf` bounds, a bound on one side only, binary constraints
  on a knapsack instance with a known optimum, binary constraints with
  a singular `H` (`nonconvex`).
- **Limits**: `time_limit` stops a branch and bound and returns
  `feasible` with an incumbent or `time_limit` without; `iter_limit`
  returns `iteration_limit`.
- **Mapping tests**: duals' signs and order, active set, statuses for
  every DAQP exit flag DAQP can return on purpose, and the `unknown`
  fallback.
- Tests self-sufficient, under `shuffle = TRUE`, CRAN suite under 15 s.

## 15. Performance targets

Measured by `tools/run-benchmarks` against `quadprog`, `osqp` and
`ROI.plugin.qpoases`, with DAQP compiled exactly as the package compiles
it (not with CMake's `-O3`), recorded in the design:

- Single solves, `n = 10…500` dense, strictly convex: faster than
  `quadprog::solve.QP()` at every size.
- A rolling sequence of 1,000 problems differing only in `f`: at least
  5× faster through `qp_model()` than 1,000 calls of `solve.QP()`.
- Overhead of the R wrapper below 20 µs per solve for `n = 10`.

## 16. Decisions

| # | Question | Decision |
| --- | --- | --- |
| D1 | Solver | DAQP, vendored |
| D2 | Precision | double only |
| D3 | Matrix layout | transpose `A` to row-major in C, once per setup |
| D4 | Statuses | result fields; errors only for invalid input and refused setup |
| D5 | quadprog compatibility | a separate function with quadprog's shape |
| D6 | PSD `H` | DAQP's automatic regularisation, reported in `regularised` |
| D7 | Warm starts | through `qp_model()` only |
| D8 | Result shape | shared with `zuprima` (`x`, `value`, `status`, …) |
| D9 | Info function | `zuqp_info()`, the package name (R1) |
| D10 | Timer | compile DAQP with `-DPROFILING`, so `time_limit` and the time fields work |
| D11 | Bounds | `ms = n` always; `±DAQP_INF` for unbounded variables |
| D12 | Model sense | `qp.sense` aliases `work.sense`; equality status fixed per model |
| D13 | Vendoring pin | a commit that includes `361ce54`, recorded in `VERSION`; a tag when one covers it |
| D14 | One-shot path | `setup_daqp_main()` + `daqp_solve()` inline, not `daqp_quadprog()`, to read `n_prox` |
| D15 | Name | `zuqp`; "DAQP" in the title and description |

Reasons where they are not in the section cited:

- D1: DAQP is the only maintained dense active-set QP solver in plain C
  with an MIT licence and no dependencies; qpOASES is C++ and already
  reachable through ROI.
- D5: overloading `qp_solve()` with `solve.QP()`'s conventions (columns
  as constraints, `meq`) would make the native interface worse; a
  separate, exact compatibility function costs little.
- D6: the first draft had zuqp retry once with proximal regularisation
  after a failed factorisation. DAQP's default already does that, in
  the factorisation, before the first solve; a second mechanism on top
  would only differ in its threshold.
- D10: without it `time_limit` is silently ignored, which contradicts
  the one promise the package makes about limits.
- D15: `daqp` was free on CRAN too, but the family prefix is what
  `zu-family.md` asks for, and the library's name belongs in the title
  where search finds it.

## 17. Open questions

1. **ROI plugin in 0.1.0?** ROI users would reach zuqp with no code
   change. Recommended: after 0.1.0, as a separate small package.
2. **Duals' convention.** DAQP signs multipliers by active side;
   `quadprog` returns non-negative `Lagrangian`. Keep DAQP's in
   `qp_solve()` and convert only in `qp_solve_quadprog()`? Recommended:
   yes; Stage 1 records the sign DAQP uses for each side.
3. **Binary constraints' scope.** Expose only `x_i ∈ {0, 1}` (bounds),
   or general `A x ∈ {bl, bu}`? Recommended: general, as DAQP has it,
   with the bound case as the documented example.
4. **Is `optimal_inexact` a success?** DAQP's flag 4 means a constraint
   is violated by more than `primal_tol` after cycling. Counting it as
   success matches "a usable `x` came back"; counting it as failure
   matches "exactly as quadprog users expect". Recommended: success,
   with the status visible, until the Maros–Mészáros runs of Stage 5
   show how often it occurs.

## 18. Roadmap

The stages, their exit criteria and their issues live in
`.agents/roadmap.md`, which is authoritative. The proposal:

| Stage | Title | Size |
| --- | --- | --- |
| 0 | Identity; vendor DAQP at a commit; measure build, timer and NOTEs | S |
| 1 | `qp_solve()`: packing, bounds, equalities, statuses | M |
| 2 | `qp_solve_quadprog()` and the quadprog oracle | S |
| 3 | `qp_model()`, `qp_update()`, warm starts | M |
| 4 | PSD and LP; binary constraints | M |
| 5 | Maros–Mészáros subset, edge cases, benchmarks | M |
| 6 | Documentation and release 0.1.0 | S |

## 19. Acceptance criteria for v0.1.0

1. Builds everywhere from package sources alone; `R CMD check
   --as-cran` clean on every CI leg; `R_init_zuqp` the only export.
2. Agrees with `quadprog` on the random-problem oracle and returns its
   exact list shape from `qp_solve_quadprog()`.
3. Every DAQP exit flag at the pinned commit maps to a documented
   status; no solve outcome raises an error.
4. `time_limit` terminates a solve on every CI platform.
5. The warm-start target of §15 is met and recorded.
6. Every mapping row of §6 has a test.
7. `src/daqp/VERSION` matches the vendored files' hashes.

## 20. What this RFC does not decide

- Whether zuqp should grow sparse support by vendoring a second solver.
  It should not; sparse belongs to `osqp`, `piqp` and `clarabel`.
- Whether zuqp replaces `quadprog` for its reverse dependencies. That
  is their maintainers' choice; zuqp makes it a rename.

## 21. Changes in the 2026-10-13 revision

Corrections after re-reading the DAQP source at `361ce54`:

- **Timer.** DAQP's timer, and with it `time_limit`, `solve_time` and
  `setup_time`, exist only under `-DPROFILING`. The first draft assumed
  they were always present. D10 adds the flag; §11, §13, §14 and §19
  test it.
- **Settings defaults.** `iter_limit` is 10,000 (not 1,000),
  `eps_prox` is −1e-6 meaning automatic (not 0), `time_limit` is 0
  meaning none (not `Inf`). `eq_reduction`, `rel_subopt`, `abs_subopt`
  and `fval_bound` were missing.
- **Exit flags.** `FEASIBLE` (5), `SOFT_OPTIMAL` (2) and
  `OVERDETERMINED_INITIAL` (−6) were not mapped; §5 now maps all twelve
  with their codes and adds an `unknown` fallback.
- **Bounds.** "`ms = n` (or the bounded count)" was wrong: the first
  `ms` constraints bound the first `ms` variables, so `ms` is always
  `n` (D11).
- **Equalities.** Set explicitly as `ACTIVE + IMMUTABLE`, as DAQP's
  docs do, instead of relying on bound-equality detection.
- **Warm starts.** The active set survives an update of `H` or `A`
  only when the problem's sense array is the workspace's own (D12);
  that fixes a row's equality status for a model's life. Models are
  not equality-reduced automatically.
- **PSD `H`.** DAQP already regularises automatically; the zuqp-side
  retry and `notes` field are replaced by a `regularised` result
  field (D6, D14).
- **Binary constraints** need a positive definite `H`; with proximal
  iterations DAQP returns `NONCONVEX`.
- **Vendoring.** Ten C files, not eleven; `v0.10.3` lacks the
  `FEASIBLE` flag and five later commits, so the pin is a commit (D13).
  `daqp_update_ldp()` lives in `utils.h`, not `api.h`.
- **quadprog shape.** `iterations` is a length-two vector and `iact`
  is ordered by addition in `quadprog`; both differences are now
  documented.
- **Name.** Resolved to `zuqp` (D15); both names were free on CRAN.
- **Open questions** renumbered; whether `optimal_inexact` counts as
  success is now question 4.

Review prepared with assistance from generative AI.
