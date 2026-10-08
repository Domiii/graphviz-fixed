# many-parallel: root cause

Sample: `tests/many_parallel/tasks_and_containers_coarse.dot` (cluster-heavy state machine). Tests: `tests/test_many_parallel.py`.

## Issue 1: not vertical (fixed in the .dot via `newrank=true`)
- Default ranking `dot1_rank` (rank.c:509, picked in `dot_rank` rank.c:528) is cluster-recursive.
  - `collapse_cluster` (rank.c:334) ranks each cluster alone, then collapses it to a leader.
  - Inter-cluster edges become soft slack pairs in `interclust1` (class1.c:33): their direction is not enforced.
  - Nested clusters without inner edges (`cluster_tasks`, `cluster_scripts`) collapse onto one rank.
- Sample: 5 ranks, 7 nodes/rank, 10 constraint edges pointing up/flat, 136 crossings.
- `newrank=true` → `dot2_rank` (rank.c:1077) ranks globally: 14 ranks, 0 upward edges.
- Other built-ins evaluated on top of newrank, none useful: `remincross`, `mclimit`, `searchsize`, `ordering=out`, `clusterrank=global`, `ratio=compress`, `nodesep`, `ranksep`, `splines=ortho` (ortho: 65 crossings).
- Minimal repro: `tests/many_parallel/cluster_rank_collapse.dot` (default xfail-strict, newrank passes).

## Issue 2: tangles left with newrank
### A. stale `ND_order` after newrank fill-node removal (dotinit.c)
- newrank: `realFillRanks` (mincross.c:975) adds NORMAL fill nodes to clusters' empty ranks; `removeFill` (dotinit.c:242) drops them after `dot_position`.
- `remove_from_rank` (dotinit.c:216) shifted the root rank array without renumbering `ND_order` or the cluster rank windows.
- `dot_splines` indexes ranks by `ND_order`: `neighbor` (dotsplines.c:2234) returns the node itself as left neighbour and skips the real right one.
  - `maximal_bbox` (dotsplines.c:2175) boxes span over neighbours; `recover_slack` (dotsplines.c:2061) then moves label vnodes onto neighbour labels.
- Effect (stock 2.43 too): 6 overlapping edge labels, braids `harness→post / heartbeat→post`, `spawn→idle / idle→drain`.
- Fix: renumber shifted nodes; `shift_cluster_ranks` (dotinit.c:200) moves/shrinks cluster windows.
- Follow-up: `checkFlatAdjacent` (flat.c:209) ignores fill nodes; else a flat edge is routed as non-adjacent while its ends are now neighbours → overlapping end boxes, `Pshortestpath failed`.
- Repro: `tests/many_parallel/newrank_fill_order.dot`.

### B. mincross local minimum: adjacent chains crossing (mincross.c)
- `retry→queued` (back-edge chain) crossed `queued→claim` several ranks below `queued`. Removing it needs a simultaneous multi-rank swap; transpose only swaps neighbours on one rank.
- Fix: `untangle_adjacent_chains` (mincross.c:1645), called from `cleanup2` (mincross.c:856) after `rec_reset_vlists` (exact cluster windows), before FLATORDER edges are dropped.
  - `untangle_fan` (mincross.c:1612): per node and direction, walk its vnode chains rank by rank. At the first gap where two neighbouring chains cross, swap their vnode prefixes.
  - Chains with the same port, plain chain vnodes only (in=out=1, no flat edges), same cluster windows; skipped with `concentrate`.
  - The segment set is unchanged except the crossing pair is uncrossed → crossings strictly drop → terminates.

### C. spline fitter backtracks (pathplan/route.c)
- `queued→claim` straight-run mode (`straight_len`, dotsplines.c:1784) ends its first piece with a vertical tangent (dotsplines.c:1804). The fitter's first candidate (`a = 4`, route.c:227) put the end control point above the start: the curve leaves queued upward, then drops across `queued→spawn`.
- `idle→ce_idle` / `ce_idle→idle`: both pieces bulge back into idle's rank band and cross twice; their polylines don't cross.
- `splinefits` accepts any candidate inside the box polygon, even if it backtracks along an axis the route is monotone in.
- Fix: `route_dirs` (route.c:289) and `splineismonotone` (route.c:306). A piece whose route is monotone in x or y must not backtrack along that axis (Bernstein derivative test); otherwise the existing loop shrinks `a` or splits.

## Rejected / kept
- Clamping overlapping flat end boxes instead of the `checkFlatAdjacent` skip: brings back `Pshortestpath failed` on 3 fuzz graphs.
- v1 of B ran before `rec_reset_vlists` (stale cluster windows → `trouble in init_rank`, nodes inside foreign clusters) and checked all pairs (O(deg² · clusters)); replaced by the fan walk.
- Fixing the default ranking (`dot1_rank`): upstream design; `newrank` is the built-in lever.
- Known limitation: adjacent flat edges with record ports (`make_flat_adj_edges`) still fail, as in base; A's skip makes such pairs adjacent more often (4/2000 fuzz graphs vs 16 fixed).
