"""
Corner compass ports (ne, nw, se, sw) under a rotated rankdir.

A port at a node corner lies on two sides. With rankdir BT/LR/RL its sides
are rotated into dot's internal top-to-bottom frame like a one-sided port's
(shapes.c invflip_side); the path then leaves the node across one of them.
"""

import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.append(os.path.dirname(__file__))
from gvtest import which  # pylint: disable=wrong-import-position
from test_edge_shapes import draw_points  # pylint: disable=wrong-import-position
from test_many_parallel import crossings  # pylint: disable=wrong-import-position
from test_port_fan import _GAP, gap, sibling_pairs  # pylint: disable=wrong-import-position

TESTS = Path(__file__).parent
RANKDIRS = ("TB", "LR", "BT", "RL")
OUTWARD = {"n": (0, 1), "s": (0, -1), "e": (1, 0), "w": (-1, 0)}  # y up

def run_dot(source: str, rankdir: str) -> dict:
    """lay out source (rankdir appended: last assignment wins) as JSON"""
    source = source[: source.rindex("}")] if "}" in source else source
    proc = subprocess.run(
        [which("dot"), "-Tjson"],
        input=f"{source}\nrankdir={rankdir};\n}}\n",
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
    assert proc.returncode in (0, 1), f"dot crashed: {proc.stderr[-500:]}"
    return json.loads(proc.stdout)


COMPASS = {"n": (0, 1), "ne": (1, 1), "e": (1, 0), "se": (1, -1), "s": (0, -1),
           "sw": (-1, -1), "w": (-1, 0), "nw": (-1, 1)}
# rankdir -> (the rankdir it mirrors, the axis it mirrors: 0 x, 1 y), or the
# quarter turn of top-to-bottom it is
SAME = {"LR": ("TB", None), "BT": ("TB", 1), "RL": ("LR", 0)}


def corner_graph(mirror: int | None) -> str:
    """corner and side ports of t (mirrored across axis: n <-> s for 1)"""
    name = {v: k for k, v in COMPASS.items()}

    def port(c: str) -> str:
        v = COMPASS[c]
        if mirror is not None:
            v = (-v[0], v[1]) if mirror == 0 else (v[0], -v[1])
        return name[v]

    edges = "".join(
        f"t:{port(c)} -> o{i}; i{i} -> t:{port(c)};\n"
        for i, c in enumerate(("ne", "se", "sw", "nw", "e", "n"))
    )
    return (
        'digraph corners {\n node [shape=box width=0.6 height=0.6 label=""];\n'
        f" t [width=1.5 height=1.5];\n{edges}"
    )


def layout_points(graph: dict) -> dict:
    """edge name -> drawn points"""
    names = {o["_gvid"]: o["name"] for o in graph["objects"]}
    return {
        f"{names[e['tail']]}->{names[e['head']]}": draw_points(e) for e in graph["edges"]
    }


@pytest.mark.parametrize("rankdir", list(SAME))
def test_corner_ports_turn_with_rankdir(rankdir: str):
    """
    A rankdir only turns or mirrors the drawing: with square nodes and no
    labels, LR is TB turned counterclockwise, BT is TB and RL is LR mirrored,
    each of the graph whose ports turn or mirror the same way. Corner ports
    lie on two sides; both move.
    """

    other, axis = SAME[rankdir]
    got = layout_points(run_dot(corner_graph(None), rankdir))
    ref = layout_points(
        run_dot(corner_graph(axis if axis is not None else None), other)
        if axis is not None
        else run_dot(corner_graph_turned(), other)
    )
    shift = None
    bad = []
    for name, points in got.items():
        if len(points) != len(ref[name]):
            bad.append(name)
            continue
        for p, q in zip(points, ref[name]):
            r = (q[1] * -1, q[0]) if axis is None else (
                (-q[0], q[1]) if axis == 0 else (q[0], -q[1]))
            d = (p[0] - r[0], p[1] - r[1])
            shift = shift or d
            if math.dist(d, shift) > 0.5:
                bad.append(name)
                break
    assert not bad, f"edges that do not follow rankdir={rankdir}: {bad}"


def corner_graph_turned() -> str:
    """corner_graph with each port turned clockwise a quarter (LR's frame)"""
    src = corner_graph(None)
    turn = {"n": "e", "ne": "se", "e": "s", "se": "sw", "s": "w", "sw": "nw",
            "w": "n", "nw": "ne"}
    return re.sub(r":([nsew]+)\b", lambda m: f":{turn[m[1]]}", src)


def sibling_metrics(graph: dict) -> tuple[int, int]:
    """crossing points and run-together pairs of edges sharing a node"""
    pairs = sibling_pairs(graph)
    cross = sum(len(crossings(p1, p2)) for _, p1, p2 in pairs)
    run = sum(1 for _, p1, p2 in pairs if min(gap(p1, p2), gap(p2, p1)) < _GAP)
    return cross, run


# repo graphs with corner ports, at 8d835b29e: sibling crossing points,
# run-together pairs, bbox width and height (pt). Neither count may grow; the
# bbox by at most _GROWTH: ports.gv LR grows 27% in height where its corner
# ports reach the side they lie on (crossings 10 -> 6, run-together 10 -> 7).
_GROWTH = 0.3
PINS = {
    ("graphs/sl_box.gv", "LR"): (0, 7, 2276, 132),  # 5 at 8d835b29e, see below
    ("graphs/sl_box.gv", "RL"): (5, 7, 2472, 132),
    ("graphs/sl_box.gv", "BT"): (0, 0, 275, 1521),
    ("graphs/st_box.gv", "LR"): (2, 6, 2337, 139),
    ("graphs/st_box.gv", "RL"): (3, 8, 2337, 139),
    ("graphs/st_box.gv", "BT"): (0, 18, 184, 1484),
    ("graphs/ports.gv", "LR"): (10, 10, 444, 392),
    ("graphs/ports.gv", "RL"): (11, 10, 444, 392),
    ("graphs/ports.gv", "BT"): (15, 14, 760, 229),
}
# sl_box LR: its self loops at corner ports pass the edge entering the node
# (7 pairs), as in TB with the ports turned a quarter (the same 7 pairs); at
# 8d835b29e the corners kept their TB sides there.


@pytest.mark.parametrize("graph_file,rankdir", list(PINS), ids=lambda v: str(v))
def test_corner_port_graph_no_worse(graph_file: str, rankdir: str):
    """
    Repo graphs with corner ports lay out no worse than at 8d835b29e.
    """

    graph = run_dot((TESTS / graph_file).read_text(encoding="utf-8"), rankdir)
    cross, run = sibling_metrics(graph)
    bb = [float(v) for v in graph["bb"].split(",")]
    width, height = bb[2] - bb[0], bb[3] - bb[1]
    base = PINS[(graph_file, rankdir)]
    print(f"{graph_file} {rankdir}: {(cross, run, round(width), round(height))}")
    assert cross <= base[0], f"sibling crossings {base[0]} -> {cross}"
    assert run <= base[1], f"run-together pairs {base[1]} -> {run}"
    assert width <= (1 + _GROWTH) * base[2], f"width {base[2]} -> {width:.0f}"
    assert height <= (1 + _GROWTH) * base[3], f"height {base[3]} -> {height:.0f}"


