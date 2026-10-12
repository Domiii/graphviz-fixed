# port-fan: root cause

Samples: `tests/port_fan/trace_two_spawns.dot`, `tests/port_fan/trace_fan_out.dot` (graph-intelligence agent traces: HTML-table node, each Agent row a `PORT` cell, edges `tailport="rN_k:e"` to subagent stubs; `rankdir=TB`). Probes: the other `tests/port_fan/*.dot`. Tests: `tests/test_port_fan.py` (`dev/test.sh tests/test_port_fan.py`).

## Commits
- ordering (`dominik/fix-port-order`): Cause 1 fix and `farface_left`; tests whose layout needs routing nesting `xfail` there.
- corner ports (`dominik/fix-invflip-side`): `invflip_side` (shapes.c:2607) turns each side of a corner port (ne, nw, se, sw) with rankdir BT/LR/RL; tests `tests/test_corner_ports.py`.
- nesting (`dominik/fix-port-fan`, `dominik/port-fan-simplify`, on the corner-port fix): Cause 2 fix; those `xfail`s pass.

## Model (internal TB frame; splines and polylines only, `dot_ports_wrap`)
- A port is a point on the node perimeter; an edge end leaves it toward one adjacent rank.
- Order toward one rank = perimeter order, cut at the far pole (`port_pos`, mincross.c:618).
  - below: counterclockwise from north; above: clockwise from south.
  - ports not on a node side leave through the face toward that rank; compass ports inside a table keep their own perimeter order as a tiebreak.
  - at the pole an end wraps away from the opposite fan (`farface_left`, splines.c:379).
- A strand: an end on a lateral side passes that side; an end on the far face wraps the face, then a side; an end at a corner of the near face is on that side at the corner, entry 0: an empty side stretch, begin/endpath's box (`run_of`, dotsplines.c:1761).
- Strands on one side toward one rank share stretches: the side down to the near corner, the far face to its corner.
- Order on a side (`strand_cmp`, :1840): perimeter distance from the near corner (`run_t.key`); one point: far end farther back = inside, then copies of one edge left to right; = the order mincross optimises.
- Objective per side (lexicographic): crossings only where the far ends' order disagrees with the ports'; consecutive strands ≥ sep = nodesep/2 apart on shared stretches; else the shortest path.

## Cause 1: mincross sees no crossing
- `local_cross`, `in_cross`, `out_cross` compared ports by `p.x`; all `:e` ports of a right-column cell share `p.x` (shapes.c:2735) → 0 crossings → `mincross` stops (mincross.c, `cur_cross == 0`) → heads in declaration order = crossing order on the east side.
- Fix: compare `port_pos` (cached per edge and call, `cached_pos`, mincross.c:637); seed `build_ranks` with compass ports in perimeter order, reversed as `build_ranks` flips the graph it ranks (`enqueue_neighbors`, mincross.c:1375; avoids O(n²) transposes on large fans).
- Median `VAL` stays `port.order`.

## Cause 2: strands share the corridor beside / around the node
- begin/endpath give every lateral port the same side box, every far-face port the same face box; shortest paths hug the node side and corner; splines bulge to the graph bound.
- Fix: channel tracks from the order, built as boxes for the stock `Pshortestpath` / `Proutespline` (no routing order, no routed-neighbour bounds):
  - index (`strands`, dotsplines.c:1913; once per node, cached per routing `strand_of` :1992, exact match): sort; k = tracks inside each strand on its side and on the far face; each of the `ED_count` copies of an edge is a strand (along a chain its far end is its lane); ports resolved as begin/endpath resolve them (`far_port` :1834).
  - track k: the band (k+½ .. k+1)·sep off n (`nest_end` :2109, `local_box` :2089); sep = min(nodesep/2, room beside n / tracks); only where ≥ 2 strands share the side.
  - entry: ends at one point spread over the port's cell, ½ sep inside its edges (`cell_keys` :1792; strands at adjacent cells ≥ sep apart; samehead / sametail ports have no side and meet); a stub joins the entry to the track, half way to the next entry inside.
  - far face: a stub rises to the face band (k+½ .. k+1)·sep_f, then the side track; sep_f from `FAR_FACE_SHARE` of the gap beyond the face; a path turns back at one end only (routespl.c:341 builds up - left - down chains mirrored).
  - leaving: strands going back across their own track inner first, the others outer first; each leaves its track at the least depth past the near face that keeps it sep clear of both straight paths of the one before it, never higher (`depths` :1881, `clear_depth` :1871); a band at that depth toward its far end, where the gap has room for it; past it, the rank box keeps off the far end of the strand leaving after it (± sep/2, where the far ends are in order; `below` :2189: the separation of consecutive strands, in the rank box both ends of a gap share).
  - copies of an edge (`make_regular_edge` :2447, `regular_path` :2216), numbered across the groups dot_splines sorts them into: no end nests: one route, copies translated by nodesep (as dot); else each copy routes along its own strands and its own lane, centred on the chain's virtual nodes, left to right at both ends and along the chain (`lane_box` :2205); the chain keeps its width (no slack recovered).
  - room at positioning, same index: `dot_lateral_room` (:2047) widens the same-rank gap by sep per track beyond the first (position.c:338); `nest_room` (position.c:127, `dot_nest_room` :2064) widens a rank gap until the bands of both its ranks' ends (tails above, heads below) are sep apart and from the ranks, and to 2·(face tracks + 1)·sep beyond a far face; gaps from `dot_rank_gap` (:2036).
  - routespl.c `checkpath` drops boxes under 0.01pt; `routesplines_` routes the boxes it keeps (a dropped box left a stale copy of the last one: "illegal values of prev 1 and next 1", a lost edge; `bundle3_ports` RL).
  - tolerances: sep, `FAR_FACE_SHARE`, `PORT_TIEBREAK` (mincross.c:44).

## Tests (`tests/test_port_fan.py`)
- lays out without crash: every `tests/port_fan/*.dot` (incl. `complex_weave_sample2`: rank -1 for a flat edge label), 4 rankdirs.
- clean in all rankdirs: both samples, 24 Agent rows, interior compass ports, same-rank neighbour, table chain, fan in a cluster, same-port fan, multi-edge bundle (one and two ranks long), bundles nested at both ends (chain: not LR, both ends on near faces; beside fans: not RL, five face wraps over ranksep/2), side ports to side ports one rank apart (except below); `port_sides` TB/LR/RL.
- no routing errors: all fixtures, `sliver_boxes` TB/BT, `flat_label_top[_ports]` (rank -1).
- `xfail` (strict; fewer failing pairs than at 847d58473, `test_known_cases_no_worse`):
  - `neighbour_heads`, `shared_heads`, `thin_box_neighbour`: an edge from another node enters a head of the fan or the table's far face; such pairs meet outside the nested strands.
    - `neighbour_heads` RL cross: `n1 → a6` pulls `a6` past `a5`; the far ends are in the opposite order of the ports: one crossing (the minimum). At 8d835b29e the pair touched at one Bézier point (90, 54), which `crossings` does not count.
  - `table_chain` RL: face wraps share half the rank gap; under 3pt apart.
- repo graphs (`test_repo_graph_no_worse`): triedds, html2, 1332, 241_1, b15, 1323 no worse than at 847d58473 (2765_1: not pinned); bbox ≤ 1.2 × at 8d835b29e (b15 +17% height: its gaps hold nested departures; run-together 40384 → 23629pt, crossings 313 → 307).
- cost (`test_big_fan_cost`): 1500 `:e` ports against the same fan without ports, best of 3: time ratio ≤ 1.2 × ×1.27, height ratio ≤ 1.01 × ×1.213 (8d835b29e: ×2.54, ×1.192).
- corner ports (`tests/test_corner_ports.py`): sl_box / st_box / ports in LR/RL/BT no worse than 8d835b29e, bbox ±1pt; pinned changes: sl_box LR (run-together 5 → 7, width 2276 → 2393; as TB with the ports turned), ports.gv LR (height 392 → 500, crossings 10 → 6, run-together 10 → 7).

## Rejected
- Bounding by routed neighbours' splines (8d835b29e `bound_side`, slicing, routing order `nest_key`): order-dependent; ~600 lines of tolerances.
- Turn by the far end across n's side, one turn track per inner turning strand: strands back across their own track but not n's side cross the inner strands' tracks; outward fans of shallow strands run together.
- One leaving depth per track for every strand: rank gaps grow by sep per strand where the paths are steep.
- Copies of a bundle offset along the path normal: they cross each other at the bends; translated: they coincide on horizontal stretches; routed apart but through one lane: they coincide along the chain.
- Polyline nudging inside tubes: `splinefits` rejects fits shorter than the input polyline; needs a clamp and a fallback.
- Inner-first / outer-first routing orders; lanes for every end; wrap side by far end; new median `VAL` (see 8d835b29e notes).
