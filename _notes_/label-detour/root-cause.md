# label-detour: root cause

## Symptom
Labeled fan-out edges from `n14` (`tests/label_detour/fanout_label_detour.dot`) swing far
left and back: `n14 -> n31` S-curves around its own label then over `n21`; `n14 -> n21` loops
to x=459 to reach its label. Stock 2.43 and upstream base `71e607d5e` show no swing.

## Verdict
Regression from the `done:` pass added in 1c09b48a5 (`dot_mincross`, lib/dotgen/mincross.c:401):
`place_free_long_virts` runs on the final order, with no later mincross/transpose.
- Its walk (`place_free_long_virts`, mincross.c:1467) moves a free `constraint=false` virt
  toward `spine_of_rank` (mincross.c:1351) = leftmost NORMAL with edges, and stops only on
  `left2right`. No bound on what it crosses.
- Multi-column ranks: "spine" is some node in the leftmost column, so the walk drags whole
  back-edge chains across unrelated columns, per rank (chain zigzags).
- Here: `n54/n56/n58 -> n26` virts dragged on r6-r12 (e.g. r8 `n58->n26` order
  16 -> 2, r10 15 -> 3), straight through `n14`'s fan. Weighted `ncross` 93 -> 235.
- `n14`'s fan inherits those crossings; x-positioning bends it around the dragged chains.
- Not spline routing: the pass is order-only and toggling it alone flips the result.

## Evidence (toggle builds, same fixture)
| build | n14->n31 | n14->n21 | n14->n53 | n14->n40 | all-edge backtrack |
|---|---|---|---|---|---|
| upstream `71e607d5e` | 0.0 | 0.0 | 11.6 | 153.0 | 118011 |
| + dotinit xpenalty/weight 1 only | 12.5 | 0.0 | 0.0 | 874.7 | 10994 |
| HEAD (`done:` pass on) | 198.8 | 490.1 | 817.1 | 1745.1 | 26902 |
| HEAD, `done:` pass off | 12.5 | 0.0 | 0.0 | 874.7 | 10994 |
| HEAD, in-`flat_reorder` call off | 198.8 | 490.1 | 817.1 | 1745.1 | 26902 |
- `done:` pass off breaks `left_stays_left` (why it was added).

## Fix
Bound the walk to its purpose (pull a free virt off a Note's far outer): it may only cross
NORMALs of the rank's constraining-flat (Note/Mid) component (`in_flat_component`,
mincross.c:1440; guard mincross.c:1506). Ranks without such a component → no move.
- fanout has no flat edges (no flat component), so there the fix just turns the pass off:
  identical to "`done:` pass off".
- All `edge_shapes`, `complex_weave`, `_notes_/long-edge-weave/fixtures`: layouts
  equal to HEAD (edge backtrack + bb); #2368 test passes.

## Also removed: constraint=false flats as hard order (same commit's weight change)
- 1c09b48a5 gives nonconstraint edges weight 1, so their flat edges pass the `ED_weight == 0`
  tests in `flat_search` and `constraining_flat_edge` (mincross.c) and become hard
  left-to-right constraints (upstream: weight 0 → none).
- Fix: `nonconstraint_flat` skips them there (walks `ED_to_orig` to the NORMAL original;
  order edges without one are never `agxget`).
- Repro `tests/nonconstraint_flat/flat_no_order.dot`: HEAD forces `y` left of `x` (1 crossing),
  fix 0 = upstream. `tests/1213-2.dot`: 2 -> 1 crossings (upstream 0). Other fixtures unchanged.
- Kept the in-`flat_reorder` `place_free_long_virts` call: dropping it changed cluster
  fixture `tests/2108.dot`, so it is not layout-neutral.

## Rejected
- Drop the `done:` pass: breaks `left_stays_left`.
- Crossing-guarded per-chain pass (one side per chain, transpose repair, keep iff `ncross`
  does not grow): fixes fanout, tests pass, but flip fixtures (intended R-then-L corridor
  needs one spine crossing) collapse to one side and wrap a Note (backtrack 94 -> 304).
- Global transpose after the pass: leaves 1000 (`left_stays_left`) / 125 (fanout) residual.

## n14 -> n40 (874.7): x-positioning trade-off, not fixed
- Final mincross order is identical for nonconstraint `weight` 1 vs 0 (rank dumps diff empty),
  so the swing is purely x-positioning on a fixed order.
- Order: `n14->n40` chain runs left of nonconstraint back-edge chains `n44->n31`,
  `n51->n53`, `n47->n53` on r14-r20 and crosses them at r20-r24 (`n40` is right of them).
  Crossing count is the same wherever it crosses.
- 1c09b48a5 gives those chains weight 1 (`dot_init_edge`, lib/dotgen/dotinit.c:73), times
  `virtual_weight` C_VV=4 (mincross.c `virtual_weight`; class2.c `make_chain`) = 3 x 4 vs
  `n14->n40`'s 4. Network simplex minimizes sum w|dx| → bends `n14->n40` (label pinned at
  r16 against the `n13->n49` chain).
- Knob sweep (fanout n40 / fanout total / sample2 total):
  | nonconstraint chain weight | n40 | fanout | sample2 |
  |---|---|---|---|
  | 1 x C_VV (HEAD, fix) | 874.7 | 10994 | 370 |
  | 1, no `virtual_weight` | 160.8 | 13414 | 981 |
  | 0 (upstream) | 120.3 | 21351 | drift (complex-weave group 2) |
- "No multiplier" re-drifts exactly 1c09b48a5's targets: sample2 `evaluateChecks->ToolRunner.run`
  122 -> 374, `finishStep->ToolRunner.run` 26 -> 248, `runAbstractSequence->observeMissingBug`
  37 -> 213; `long-edge-weave/clean-flats` 69 -> 121. Tests still pass, quality regresses.
- Only the smallest setting (1, no multiplier) unbends `n40` here; more chains would beat it too.
  No robust weight setting → inherent to the 1c09b48a5 weighting. Real fix would be
  order-level (where `n14->n40` crosses the bundle), out of scope.

## Visual check (after.svg)
- n14 fan: no loops/S-curves; labels clear of edges/nodes.
- Remaining, not caused by this pass's guard (same as "`done:` pass off"):
  - `n24 -> n10` swings right (482 -> 1123): its chain sits right of `n26` on r3-r5 and
    crosses `n02->n39`, `n26->n32`, `n58->n26` into `n10` (left). Mincross local order;
    HEAD's cross-column drag happened to help it.
  - `n44 -> n31` mild hook below `n31` (598 -> 742).

## Limitations
- Fanout bb widens 2744 -> 3236 (the drag had compressed it); upstream base is 3776.
- A free virt separated from the Note by an unrelated node stays on the far outer.
