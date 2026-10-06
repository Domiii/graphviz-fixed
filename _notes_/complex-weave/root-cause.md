# complex-weave: root cause

## Symptom
Labeled `constraint=false` back edges (e.g. `Low -> Top`, `evaluateChecks -> ToolRunner.run`,
`resolve_missing_email -> runEmailStep`) wrap the far side of an unrelated column, widen the
graph, and cross unrelated edge labels.

## Verdict
Root cause is in `dot_init_edge` (lib/dotgen/dotinit.c): nonconstraint edges get
`ED_xpenalty = 0` and `ED_weight = 0`.
- xpenalty 0 → mincross ignores crossings on the back edge's virtual chain, so the chain is
  parked at the rank end (far side of the other column).
- weight 0 → x-positioning does not straighten the chain, so it meanders (drift right).

The prior fix (`place_free_long_virts`, 2c2f9bdc0) is a post-hoc side correction and cannot fix
this group.

## Bug groups
1. **Far-side parking (ordering).** `far_side_back_edge`, `sample2_min` evaluateChecks /
   runRequirementsJudge / resolve_missing_email. Final order was `L-col, R-col, V[Low->Top]` on
   every rank. Cause: xpenalty 0.
2. **Chain drift (positioning).** With xpenalty 1 but weight 0 (`EXP2`),
   `V[runRequirementsJudge->ToolRunner.run]` x goes 207→241→268 and still exceeds the band.
   Cause: weight 0.
3. **Prior-fix interaction.** Once xpenalty > 0, `is_free_long_edge_virt` (xpenalty==0 test)
   rejects the virt, and later mincross passes undo `place_free_long_virts`
   (`left_stays_left` breaks: `Mid Note V` instead of `V Mid Note`).

## Why the prior fix missed (evidence, GV_WEAVE_DEBUG dumps)
- `flat_reorder` returns early: `FLAT_REORDER g=TaskJourney has_flat=0`,
  `g=FarSideBackEdge has_flat=0`. So `place_free_long_virts` never ran.
- Forcing it (`EXP1`/`EXP7`) still does nothing: every candidate is `accept=1` but `side=0`
  (`place V[Low->Top] side=0 spine_ord=0 before=2 after=2`). Labels double ranks; odd ranks are
  virtual-only (`PFLV r1 spine=NULL` … r13), so neighbor side detection finds no spine, and
  there are no flats for the fallback. `resolve_*` label virts live only on no-spine r11.
- Widening `is_free_long_edge_virt` (`EXP5`) is irrelevant for these fixtures (all already
  `xpen=0 accept=1`).

## Experiment matrix (pytest: complex_weave 6 xfails + edge_shapes 5 + #2368 2)
| env | complex_weave fixed | edge_shapes | #2368 |
|---|---|---|---|
| EXP1 / EXP5 / EXP7 / combos | 0/6 | 5/5 | ok |
| EXP2 (xpenalty 1) | 5/6 | left_stays_left FAIL | ok |
| EXP4 (weight kept) | 1/6 | 5/5 | ok |
| EXP6 (xpenalty 1, weight 1) | 6/6 | left_stays_left FAIL | ok |
| EXP2+EXP5+EXP7 | 5/6 | 5/5 | ok |
| **EXP6+EXP5+EXP7** | **6/6** | **5/5** | ok |

## Ranked solutions
1. **(v1, shipped)** nonconstraint edges keep xpenalty 1 / weight 1; free-virt test uses
   `nonconstraint_edge(orig)`; rerun `place_free_long_virts` after all mincross passes.
   Small, evidenced, fixes all groups. Risk: global layout change for every
   `constraint=false` edge.
2. Scope (1) to back edges only (tail rank > head rank) or labeled ones. Less churn, but
   `dot_init_edge` runs before ranking, so needs moving to class2 / after rank.
3. Make `place_free_long_virts` corridor-aware (spine per column, walk chain across
   no-spine ranks). Fixes ordering without touching penalties but not group 2 drift;
   larger and more heuristic.
4. Low weight (e.g. weight 1 but keep xpenalty 0) — insufficient (EXP4: 1/6).

## BEST v1
Solution 1. Files: `lib/dotgen/dotinit.c` (`dot_init_edge`), `lib/dotgen/mincross.c`
(`orig_edge_of_virt`, `is_free_long_edge_virt`, `dot_mincross` `done:` pass).

## v1 limitations
- Changes layout of every graph with `constraint=false` edges (reference-output churn).
- User `weight=` on `constraint=false` edges is overridden to 1 (upstream forced 0).
- Final `place_free_long_virts` runs on the root after cluster mincross; `side` still
  derives from the leftmost-NORMAL "spine", not columns.
- `place_free_long_virts` itself is a no-op on these fixtures; the fix is the penalty/weight.
