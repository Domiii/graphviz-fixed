"""
Edges fanning out of ports along one node side.

An HTML-table node whose rows each carry a port on the same side (`tailport=
"row:e"`) fans edges out to nodes on the next rank. Such ports share one x, so
port x alone cannot order their edges; on the east side declaration order is the
crossing order (the top port's edge is outermost). The edges share the corridor
beside the node and nest there; a port facing away from the other end wraps
around the node, and those wraps nest too.
"""

# Do not overfit: assert only the obvious, only-correct behavior (two edges sharing
# a node, attached at distinct ports, neither cross nor run together; layout
# without routing errors). Repo graphs with crossings of their own must not get
# worse than at 847d58473.

import json
import math
import os
import re
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import pytest

sys.path.append(os.path.dirname(__file__))
from gvtest import which  # pylint: disable=wrong-import-position
from test_edge_shapes import draw_points  # pylint: disable=wrong-import-position
from test_many_parallel import (  # pylint: disable=wrong-import-position
    bezier_polyline,
    crossings,
)

TESTS = Path(__file__).parent
FIXTURES = TESTS / "port_fan"

_GAP = 3.0  # closer than this, two strokes read as one
_ERRORS = re.compile(
    r"triangulation failed|Pshortestpath|not in any triangle|"
    r"cannot find triangle path|lost .* edge|illegal values|Error"
)
RANKDIRS = ("TB", "LR", "BT", "RL")

# fixture -> rankdirs whose layout can be crossing-free with dot's routing
CLEAN = {
    "trace_two_spawns.dot": RANKDIRS,  # graph-intelligence sample: 2 spawns
    "trace_fan_out.dot": RANKDIRS,  # sample: 5 spawns over 3 responses
    "agent_rows_24.dot": RANKDIRS,  # many Agent rows
    "port_sides.dot": ("TB", "LR", "RL"),  # in/out ports on every side
    "interior_ports.dot": RANKDIRS,  # compass ports inside the table
    "same_rank_neighbour.dot": RANKDIRS,  # a node beside the fan
    "neighbour_heads.dot": RANKDIRS,  # that node feeds heads of the fan
    "shared_heads.dot": RANKDIRS,  # several nodes beside it feed heads
    "thin_box_neighbour.dot": RANKDIRS,  # a thin box beside the fan
    "table_chain.dot": RANKDIRS,  # a table below a head fans out too
    "fan_in_cluster.dot": RANKDIRS,  # the fan inside a cluster
}
# fixture -> rankdirs that must lay out without routing errors
ROUTES = {
    **{f: RANKDIRS for f in CLEAN},
    "sliver_boxes.dot": ("TB", "BT"),  # boxes under 0.5pt tall
    "flat_label_top.dot": RANKDIRS,  # rank -1 for a flat edge label
    "flat_label_top_ports.dot": RANKDIRS,  # ... and a port fan
}

# Cases that still fail, with sibling pairs that cross / run together as
# measured at base commit 847d58473; neither count may grow.
_FOREIGN = (
    "an edge from another node enters a head of the fan, or reaches the "
    "table's far face where the fan wraps: such pairs meet outside the "
    "nested port runs"
)
_TIGHT = "wraps over the far face share half the rank gap: under 3pt apart"
KNOWN = {
    ("neighbour_heads.dot", "TB"): ((6, 9), _FOREIGN, ("cross", "run")),
    ("neighbour_heads.dot", "LR"): ((0, 1), _FOREIGN, ("run",)),
    ("neighbour_heads.dot", "BT"): ((5, 7), _FOREIGN, ("cross", "run")),
    ("neighbour_heads.dot", "RL"): ((7, 8), _FOREIGN, ("run",)),
    ("shared_heads.dot", "RL"): ((5, 6), _FOREIGN, ("cross", "run")),
    ("thin_box_neighbour.dot", "RL"): ((3, 3), _FOREIGN, ("cross", "run")),
    ("table_chain.dot", "RL"): ((6, 7), _TIGHT, ("run",)),
    ("port_sides.dot", "RL"): ((10, 13), None, ()),
    # in-edge ports below out-edge ports: dot routes no crossing-free layout
    ("port_sides.dot", "BT"): ((19, 17), None, ()),
}

# repo graphs this change affects, as measured at base commit 847d58473:
# sibling crossings, run-together length (pt), routing error lines; none may
# grow
REPO = {
    "../graphs/directed/triedds.gv": (2, 45, 0),
    "graphs/html2.gv": (0, 0, 0),
    "1332.dot": (29, 302, 0),
    "241_1.dot": (5, 263, 0),
    "graphs/b15.gv": (322, 40588, 0),
    "1323.dot": (37, 3256, 2),  # lost flat edges between record ports
}

# every fixture lays out in each rankdir (complex_weave_sample2: no ports, a
# flat edge label on the top rank)
LAYOUTS = sorted(FIXTURES.glob("*.dot"))


def run_dot(path: Path, rankdir: str | None = None) -> tuple[dict, str]:
    """lay out path (rankdir: last assignment wins) as JSON; routing errors"""
    source = path.read_text(encoding="utf-8", errors="replace")
    if rankdir:
        source = source[: source.rindex("}")] + f"rankdir={rankdir};\n}}\n"
    proc = subprocess.run(
        [which("dot"), "-Tjson"],
        input=source,
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
    assert proc.returncode in (0, 1), f"dot crashed: {proc.stderr[-500:]}"
    errors = [line for line in proc.stderr.splitlines() if _ERRORS.search(line)]
    return json.loads(proc.stdout), errors


def point_segment(p: tuple, a: tuple, b: tuple) -> float:
    vx, vy = b[0] - a[0], b[1] - a[1]
    length = vx * vx + vy * vy
    t = 0.0 if length == 0 else ((p[0] - a[0]) * vx + (p[1] - a[1]) * vy) / length
    t = max(0.0, min(1.0, t))
    return math.hypot(p[0] - a[0] - t * vx, p[1] - a[1] - t * vy)


def shared_ends(p: list, q: list) -> list:
    return [x for x in (p[0], p[-1]) if any(math.dist(x, y) < 1 for y in (q[0], q[-1]))]


def gap(p: list, q: list) -> float:
    """closest approach of polyline p to polyline q, away from shared end points"""
    ends = shared_ends(p, q)
    return min(
        (
            min(point_segment(x, a, b) for a, b in zip(q, q[1:]))
            for x in p
            if all(math.dist(x, e) > 4 * _GAP for e in ends)
        ),
        default=math.inf,
    )


def run_length(p: list, q: list) -> float:
    """length of p within _GAP of q, away from shared end points"""
    ends = shared_ends(p, q)
    total = 0.0
    for a, b in zip(p, p[1:]):
        m = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        if all(math.dist(m, e) > 4 * _GAP for e in ends) and min(
            point_segment(m, c, d) for c, d in zip(q, q[1:])
        ) < _GAP:
            total += math.dist(a, b)
    return total


def sibling_pairs(graph: dict) -> list:
    """pairs of edges sharing a node, as (names, polyline, polyline)"""
    names = {o["_gvid"]: o["name"] for o in graph.get("objects", [])}
    edges = []
    for e in graph.get("edges", []):
        if "pos" not in e:  # lost
            continue
        name = f"{names[e['tail']]}→{names[e['head']]}"
        polyline = bezier_polyline(draw_points(e))
        edges.append((name, {e["tail"], e["head"]}, polyline))
    return [
        (f"{n1} / {n2}", p1, p2)
        for (n1, ends1, p1), (n2, ends2, p2) in combinations(edges, 2)
        if ends1 & ends2
    ]


def cases(table, check: str | None = None) -> list:
    """(fixture, rankdir) params; xfail where KNOWN says check still fails"""
    out = []
    for f, rs in table.items() if isinstance(table, dict) else table:
        for r in rs if isinstance(table, dict) else (rs,):
            known = KNOWN.get((f, r))
            marks = ()
            if known and check in known[2]:
                marks = pytest.mark.xfail(strict=True, reason=known[1])
            out.append(pytest.param(f, r, id=f"{f}-{r}", marks=marks))
    return out


@pytest.mark.parametrize(
    "path,rankdir",
    [
        pytest.param(p, r, id=f"{p.name}-{r}")
        for p in LAYOUTS
        for r in RANKDIRS
    ],
)
def test_lays_out(path: Path, rankdir: str):
    """
    dot lays the graph out without crashing.
    """

    run_dot(path, rankdir)


@pytest.mark.parametrize("fixture,rankdir", cases(ROUTES))
def test_routes_without_errors(fixture: str, rankdir: str):
    """
    Every edge is routed: no path planning errors, no lost edges.
    """

    _, errors = run_dot(FIXTURES / fixture, rankdir)
    assert not errors, f"routing errors: {errors}"


@pytest.mark.parametrize("fixture,rankdir", cases(CLEAN, "cross"))
def test_sibling_edges_do_not_cross(fixture: str, rankdir: str):
    """
    Two edges sharing a node do not cross.
    """

    graph, _ = run_dot(FIXTURES / fixture, rankdir)
    bad = [name for name, p1, p2 in sibling_pairs(graph) if crossings(p1, p2)]
    assert not bad, f"sibling edges cross: {bad}"


@pytest.mark.parametrize("fixture,rankdir", cases(CLEAN, "run"))
def test_sibling_edges_do_not_run_together(fixture: str, rankdir: str):
    """
    Two edges sharing a node keep apart away from shared end points.
    """

    graph, _ = run_dot(FIXTURES / fixture, rankdir)
    bad = [
        f"{name} ({d:.1f}pt)"
        for name, p1, p2 in sibling_pairs(graph)
        if (d := min(gap(p1, p2), gap(p2, p1))) < _GAP
    ]
    assert not bad, f"sibling edges run together: {bad}"


@pytest.mark.parametrize("fixture,rankdir", cases(list(KNOWN), "known"))
def test_known_cases_no_worse(fixture: str, rankdir: str):
    """
    Where sibling edges still cross or run together, fewer pairs do than at
    847d58473.
    """

    graph, _ = run_dot(FIXTURES / fixture, rankdir)
    pairs = sibling_pairs(graph)
    cross = sum(1 for _, p1, p2 in pairs if crossings(p1, p2))
    run = sum(1 for _, p1, p2 in pairs if min(gap(p1, p2), gap(p2, p1)) < _GAP)
    base_cross, base_run = KNOWN[(fixture, rankdir)][0]
    assert cross <= base_cross, f"crossing sibling pairs {base_cross} -> {cross}"
    assert run <= base_run, f"run-together sibling pairs {base_run} -> {run}"


@pytest.mark.parametrize("graph_file", REPO, ids=lambda f: Path(f).name)
def test_repo_graph_no_worse(graph_file: str):
    """
    Repo graphs with ports lay out with no more routing errors, sibling edge
    crossings and run-together than at 847d58473.
    """

    graph, errors = run_dot(TESTS / graph_file)
    base_cross, base_run, base_errors = REPO[graph_file]
    assert len(errors) <= base_errors, f"routing errors: {errors}"
    pairs = sibling_pairs(graph)
    cross = sum(len(crossings(p1, p2)) for _, p1, p2 in pairs)
    run = sum(run_length(p1, p2) for _, p1, p2 in pairs)
    assert cross <= base_cross, f"sibling crossings {base_cross} -> {cross}"
    assert run <= base_run + _GAP, f"sibling run-together {base_run} -> {run:.0f}pt"
