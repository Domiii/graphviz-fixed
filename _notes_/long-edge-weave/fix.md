# Fix: keep free long-edge VIRTUAL off Notes’ far outer

## HEAD on these fixtures

- These fixtures have flat edges, so the `!GD_has_flat_edges` return does not apply
- #2368’s `valid_rank` bail-out runs when a node has a constraining flat in-edge, including a head already placed by `postorder`, and the rewrite is skipped
- Notes often land on the wrong side of Mid
- The `constraint=false` edge runs straight through the Mid|Note gap
- That is not the 2.43 far-outer weave

## Far-outer dump is conditional

- True when `flat_reorder` actually rewrites order: the free long-edge `VIRTUAL` can land on the Note’s far outer
- Not the explanation of current HEAD on these fixtures
- Rewrite-case witness: [`fixtures/clean-flats.dot`](fixtures/clean-flats.dot), [`fixtures/clean-flats-polyline.plain`](fixtures/clean-flats-polyline.plain)

## Rejected

- HoldAndReinsert stays deleted. Restoring it reintroduces [Graphviz #2368](https://gitlab.com/graphviz/graphviz/-/issues/2368) (routesplines illegal values)
- Slot-fill plus `ED_xpenalty` on the virtual chain is blocked by #2368 and is not in the tree

## Rewrite / 2.43 case, not HEAD: Bound-side (“outside all NORMALs on free Bound”)

- Directionally often right (don’t park past Notes) — **not** a universal placement law
- Fails when free Bound flips per rank, or when another NORMAL already occupies far Bound

### Flip: [`fixtures/flip-sides.dot`](fixtures/flip-sides.dot) ([`.plain`](fixtures/flip-sides.plain))

- MidA: `NoteL` left of Spine; MidB: `NoteR` right
- Back-edge: `xmin=0` / `xmax=4.46`, span `4.46`, `revs=2`
- Per-rank free Bound is right then left → Bound-side still flips → full-width `rank_box` weave
- Need **chain-stable** side-of-Spine, not per-rank Bound

### Far Bound worse than gutter: [`fixtures/both-sides.dot`](fixtures/both-sides.dot)

- Rank: `LeftNote | Mid | RightBlock`
- Rewrite-case path near Mid: `x=0` (Notes far outer) — same dump
- Bound-side past `RightBlock` detours; better: gutter `Mid | RightBlock`
- Free **direction** ≠ far **Bound**

### Same: [`fixtures/far-bound-bad.dot`](fixtures/far-bound-bad.dot)

- `Note | Mid | RightKeep`; rewrite-case `x=0`; Bound past `RightKeep` wraps keep-out
- Near gutter `Mid | RightKeep` is the free-side slot

### Minimal dump: [`fixtures/one-flat.dot`](fixtures/one-flat.dot)

- `Note | Mid`; path `x=0` — smallest Notes-outer dump

## Secondary / out of scope

- `maximal_bbox` widen for long-edge dummies — workaround; does not undo order dump
- Length-oracle side pick — larger design; not required by these fixtures
