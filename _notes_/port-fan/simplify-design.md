# port-fan: one nesting mechanism (design, revised after review)

Scope: simplifies the routing nesting and room reservation of 8d835b29e (dotsplines.c +755, position.c +61, routespl.c +32). Keeps the ordering commit 450b81901 (`port_pos`, `enqueue_neighbors`, `farface_left`).

## Problem model
- Strand: one edge's path between its end at node n and the adjacent rank (tail: below; head: above).
- Obstacle: n's box. Siblings at n toward one rank pass the same side of n; their shortest paths hug the same boundary and corners → crossings at ties, run-together.
- Shared stretches (geometry):
  - side: from a strand's port along n's side to the near corner; shared with every sibling whose port lies between it and that corner.
  - far face: from a far-face port (tail TOP, head BOTTOM) along that face to the corner of its side (`farface_left`).
  - turn: a strand turns at the near corner iff its straight path from that corner, n inflated by the offsets of the strands inside it, to its far end cuts n's box; it then runs along n's near face.
- Objective, per node and side (lexicographic):
  1. no crossing between strands whose order agrees at both ends of their shared stretch; one where it disagrees (minimum, metro-line crossing minimization with fixed terminal order).
  2. consecutive strands on a shared stretch ≥ sep apart.
  3. otherwise minimal deviation from the shortest path.

## Established formulations
- VLSI channel routing, left-edge algorithm (Hashimoto & Stevens 1971): tracks = channel density.
- Metro-line crossing minimization with fixed terminals (Benkert et al. 2006; Nöllenburg 2009; Pupyrev, Nachmanson, Bereg, Holroyd 2011): the order of two paths on a shared subpath is their order where they split; minimum crossings.
- Ordered bundles in layered layouts (Pupyrev, Nachmanson, Kaufmann 2010; MSAGL).
- Nudging with separation constraints (Wybrow, Marriott, Stuckey 2009; VPSC, Dwyer, Marriott, Stuckey 2006): min Σ (xᵢ − xᵢ*)² s.t. xᵢ₊₁ − xᵢ ≥ sep. With all xᵢ* on the boundary and kᵢ strands inside i on the same stretch, the optimum is xᵢ = kᵢ·sep.
- Virtual-node waypoints (dot's own mechanism): see "Alternatives".

## Invariant (ordering ↔ routing)
- Siblings at n on one side sort by `port_pos` toward their rank, ties by far-end order (`ND_order`): the order mincross optimises (450b81901) is the order routing nests.

## Mechanism: derived channel tracks in the stock pipeline
- Tracks are not hand-placed lanes (graphviz-fixed.md): each follows from the order invariant and the shared-stretch geometry; nothing depends on routing order or on other edges' splines.
- Per node, once before routing, O(deg log deg):
  - sort siblings by the invariant; weight each by `ED_count` (a bundle of c copies takes c tracks).
  - per strand and stretch: k = weight of the strands inside it on that stretch (side: ports nearer the corner; turn: inner strands that turn; face: face strands nearer the corner), m = total on that stretch.
  - same port: strands at one port take entry points spread along the port's cell edge, sep apart in order, sep shrunk to the cell (like samehead): distinct entries, so their tracks are disjoint.
- Track of strand k on a stretch: the band [(k+½)·sep, (k+1)·sep] off n (width sep/2, sep/2 to the next track), as the corridor beside that stretch: n inflated by (k+½)·sep there (inside bound), capped at (k+1)·sep (outside bound); applied only when m ≥ 2, so a lone strand routes as in base.
  - side: side box [side + (k+½)sep, side + (k+1)sep] down to the near corner (cap only along the shared part of the side).
  - turn: rank-box band below n [nbot − (k+1)sep, nbot − (k+½)sep] up to the strand's departure x, then the descent box toward its head, capped at the outer strand's head (−sep/2) where the two heads are in order (order agrees); where it disagrees no cap: one crossing (objective 1).
  - far face: a stub from the port to the face band [face + (k+½)sep, face + (k+1)sep], then the side track.
  - head ends: mirrored (strands arriving from above).
- sep = min(nodesep/2, room/m) per stretch; room from the corridor bound beside n, the rank gap below/above n.
- `Pshortestpath` returns the track path (the VPSC optimum on the stretch); `Proutespline` fits inside the corridor unchanged: no ε test, no barriers, no clamp, no fallback.
- Bundles: `make_regular_edge` spreads copies by the track pitch (sep), within their c tracks, not by Multisep.
- Room (positioning): the same counts drive `dot_lateral_room` (beside n: m·sep) and `nest_room` (rank gap: turn tracks m·sep; a far face gets its share of the gap, so 2·m·sep).
- routespl.c: box chains up - left - down (far-face stubs going left; head-end stubs) are built mirrored (kept from 8d835b29e: the polygon builder assumes up - right - down).

## Alternatives
- Polyline nudging + tube barriers (first draft): `splinefits` (route.c:246) rejects fits shorter than the input polyline, which it takes as the shortest path; a nudged polyline isn't that inside its tube → split pieces, kinks; two-point pieces skip the barrier check (route.c:263-273); needs a clamp and a fallback that drops the guarantee, and a routespl.c API change. Rejected.
- Virtual-node waypoints: most principled (dot routes through vnodes already), but needs an intermediate rank (rank doubling changes every port layout) and separates strands only at that rank, not at n's corners. Not now; its order index is the one above.
- Bounding by routed neighbours' splines (8d835b29e `bound_side`): needs a routing order (`nest_key`) and slicing; replaced by tracks from the order.

## The 12 xfail cases
- `table_chain` RL, `port_sides` RL (_TIGHT, 2.1-2.2pt): no claim. They pass only if the face rule's 2·m·sep gap holds and sep stays ≥ 6pt; else they stay xfail with their pins.
- `neighbour_heads`, `shared_heads`, `thin_box_neighbour` (_FOREIGN, 10): out of scope (head-end spreading changes every graph with converging in-edges: own design and PR). xfails and pins unchanged.

## Covered by 8d835b29e, and here
- clearance to routed neighbours anywhere in the tail-rank boxes (stub, side, filler, first rank box) → tracks on n's stretches + the descent cap toward the head; beyond those (a neighbour's corner, vnode chains) plain shortest paths, as both.
- inner-bounding when the two ends' orders conflict → no cap there, one crossing.
- lanes only when squeezed → tracks everywhere a stretch is shared (m ≥ 2).
- routing order by nest class → dropped (corridors are deterministic).
- may regress: open-space fans where bulging fits looked smoother; strands whose heads sit between two ports' tracks. The benchmark decides.

## Code size
- Removed (dotsplines.c): `bound_side`, `bound_by`, `slice_box`, `clear_turn`, `clear_interval`, `routed_polyline`, `simplify`, `sibling_t`, `nest_key` + edgecmp key, `NEST_SHIFT`, `NEST_BIAS`, `compass_end`, `far_gap` (room from rank coordinates in one place), and constants `NEST_SAMPLES`, `NEST_SIMPLIFY`, `NEST_SLIVER`, `NEST_CONTACT`, `NEST_MIN_BOX`, `NEST_OVERLAP`, `NEST_MERGE`, `NEST_END_ON_EDGE`, `NEST_INSET`, `NEST_FACE_LEAD`, `MAX_SLICES`.
- Rewritten: `nest_end` → track corridors from the precomputed index (`run_of`, `turns` geometric, `run_count`, `nests` → one sort).
- Kept: `dot_ports_wrap`, `dot_lateral_room`, `nest_room`, `dot_nest_sep`, `below`, routespl.c mirroring, `FAR_FACE_SHARE`.
- Estimate: ~300-350 lines vs ~850; tolerances: sep (and `FAR_FACE_SHARE`).

## Complexity
- Order: O(deg log deg) per node once; corridor per edge O(1) boxes.
- mincross: cache `port_pos` per edge per `in_cross`/`out_cross` call (atan2 per pair otherwise).
- Measured, not estimated: a pytest case for the 1500-row fan (time, bounding box height).

## Risks
- Track corridors and head positions: a descent cap from an outer head that is not in order would close the corridor → cap only where orders agree.
- Same-port spreading moves endpoints inside the port's cell (visible on record/HTML ports).
- (k+½)·sep offsets keep arrows and strokes off n (m ≥ 2 only).
- Taller layouts from `nest_room`; measured in the perf case.

## Benchmark and acceptance
- `tests/test_port_fan.py` (all fixtures × 4 rankdirs): routing errors, crossings, run-together (< 3pt), no-worse pins at 847d58473, repo graphs (triedds, html2, 1332, 241_1, b15, 1323); new fixtures: same-port fan, multi-edge bundle; perf case (1500-row fan).
- Earlier fix suites (complex_weave, edge_shapes, label_detour, many_parallel, nonconstraint_flat); regression subset 241/1323/1332/2765/b15/2416.
- ASan pass on the round-2 crash probes (/tmp/pfrev/r2, /tmp/pfrev/cp) at the end (needs an ASan build dir; `dev/build.sh` has no sanitizer switch).
- Accept: every passing case passes; xfails ≤ 12, each pinned; code smaller; 1500-row fan within base + 20%.
