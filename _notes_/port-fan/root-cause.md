# port-fan: root cause

Samples: `tests/port_fan/trace_two_spawns.dot`, `tests/port_fan/trace_fan_out.dot` (graph-intelligence agent traces: HTML-table node, each Agent row a `PORT` cell, edges `tailport="rN_k:e"` to subagent stubs; `rankdir=TB`). Probes: the other `tests/port_fan/*.dot`. Tests: `tests/test_port_fan.py` (`dev/test.sh tests/test_port_fan.py`).

## Commits
- ordering (`dominik/fix-port-order`): Cause 1 fix and `farface_left`; tests whose layout needs routing nesting `xfail` there.
- nesting (`dominik/fix-port-fan`): Cause 2 fix; those `xfail`s pass.

## Model (internal TB frame; splines and polylines only, `dot_ports_wrap`)
- A port is a point on the node perimeter; an edge end leaves it toward one adjacent rank.
- Order toward one rank = perimeter order, cut at the far pole (`port_pos`, mincross.c:618).
  - below: counterclockwise from north; above: clockwise from south.
  - ports not on a node side leave through the face toward that rank; compass ports inside a table keep their own perimeter order as a tiebreak.
  - at the pole an end wraps away from the opposite fan (`farface_left`, splines.c:379).
- A run: an end on a lateral side passes that side; an end on the far face wraps the face, then a side (`run_of`, dotsplines.c:1782).
- Runs on one side toward one rank nest: nearer the corner where they leave = inside; one port: far end farther back = inside (`nests`, :1806).
- A run turns if its far end lies back across the side (`turns`, :1813).

## Cause 1: mincross sees no crossing
- `local_cross`, `in_cross`, `out_cross` compared ports by `p.x`; all `:e` ports of a right-column cell share `p.x` (shapes.c:2735) → 0 crossings → `mincross` stops (mincross.c, `cur_cross == 0`) → heads in declaration order = crossing order on the east side.
- Fix: compare `port_pos`; seed `build_ranks` with compass ports in perimeter order, reversed as `build_ranks` flips the graph it ranks (`enqueue_neighbors`, mincross.c:1357; avoids O(n²) transposes on large fans).
- Median `VAL` stays `port.order`.

## Cause 2: runs share the corridor beside / around the node
- begin/endpath give every lateral port the same side box, every far-face port the same face box; shortest paths hug the node side and corner; splines bulge to the graph bound.
- Fix:
  - room at positioning: `dot_lateral_room` (dotsplines.c:1909) widens the same-rank gap by one nest separation (nodesep/2) per sibling beyond the first (position.c:327); `nest_room` (position.c:125, from `dot_nest_runs`) widens rank gaps where turning or face-wrapping siblings nest; ranks indexed from `GD_minrank` (rank -1 of flat edge labels).
  - routing order (`nest_key`, :1843; edgecmp :597): turning runs outer first, others inner first.
  - corridor (`nest_end`, :1974): reach out no farther than the far end needs; one separation per inner sibling beside the node, beyond it for turning ones, over the face for face-wrapping ones (the face's share of the real gap, `far_gap` :2412); a stub joins the port; squeezed beside a neighbour, a compass end keeps to its own separation; a path turns back at one end only.
  - clearance (`bound_by` / `bound_side` :2342, `slice_box` :2266, `clear_turn` :2229, `clear_interval` :2122): slice the end boxes where the nearest routed nesting neighbours bend (`routed_polyline` :2192); keep the clear interval that continues the path from its start; at an end turning inside a box, keep the clear rectangle beside the end.
  - tolerances named once, one per meaning: `NEST_*`, `FAR_FACE_SHARE`, `MAX_SLICES` (dotsplines.c:51ff), `PORT_TIEBREAK` (mincross.c:44).
  - routespl.c:341: box chains that go up - left - down are built mirrored.

## Tests (`tests/test_port_fan.py`)
- lays out without crash: every `tests/port_fan/*.dot` (incl. `complex_weave_sample2`: rank -1 for a flat edge label), 4 rankdirs.
- clean in all rankdirs: both samples, 24 Agent rows, interior compass ports, same-rank neighbour, table chain, fan in a cluster (except below); `port_sides` TB/LR/RL.
- no routing errors: all fixtures, `sliver_boxes` TB/BT, `flat_label_top[_ports]` (rank -1).
- `xfail` (strict; fewer failing pairs than at 847d58473, `test_known_cases_no_worse`):
  - `neighbour_heads`, `shared_heads`, `thin_box_neighbour`: an edge from another node enters a head of the fan or the table's far face; such pairs meet outside the nested port runs.
  - `table_chain` RL, `port_sides` RL: face wraps share half the rank gap; 2.1-2.2pt apart.
- repo graphs (`test_repo_graph_no_worse`): triedds, html2, 1332, 241_1, b15, 1323 no worse than at 847d58473 (2765_1 lays out by the routing order only: not pinned).

## Rejected
- Inner-first only: inner splines sag under the node where outer ones must pass.
- Outer-first only: non-turning outer runs pass low over inner heads.
- Bounding by a one-sided cap (left / right of the sibling): wrong inside U-shaped face wraps; replaced by the clear interval.
- Merging slices by running min / to a fixed count without keeping contact: closes corridors.
- Wrap side by far end: disagrees with mincross's perimeter cut.
- Lanes for every end: tight splines run together in open space.
- New median `VAL`: triedds sibling run 45 → 127pt.
