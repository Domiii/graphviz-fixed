"""
Cluster-heavy state machine tests (many parallel same-rank nodes).

With the default ranking (newrank=false) every cluster is ranked on its own and
inter-cluster edges are soft. A nested cluster whose nodes share no inner edge
collapses onto a single rank, so a request/response chain ping-ponging between
clusters folds into a few overfull ranks: constraint edges point up or run flat,
and edges into the same node braid. newrank=true ranks the graph globally; the
sample uses it.
"""

# Do not overfit: assert only the obvious, only-correct behavior (in a TB layout a
# constraint edge points down; two edges sharing a node never need to cross; edge
# labels never overlap).

import json
import os
import re
import sys
from itertools import combinations
from pathlib import Path

import pytest

sys.path.append(os.path.dirname(__file__))
from gvtest import dot  # pylint: disable=wrong-import-position
from test_edge_shapes import (  # pylint: disable=wrong-import-position
    bbox,
    center_y,
    draw_points,
)

FIXTURES = Path(__file__).parent / "many_parallel"

_SLACK = 18.0

SAMPLE = "tasks_and_containers_coarse.dot"
COLLAPSE = "cluster_rank_collapse.dot"


def load_layout(fixture: str, newrank: bool = False) -> dict:
    path = FIXTURES / fixture
    assert path.exists(), f"unexpectedly missing test case {path}"
    source = path.read_text(encoding="utf-8")
    if newrank:
        source = re.sub(r"\{", "{ newrank=true;", source, count=1)
    return json.loads(dot("json", source=source))


def nodes(graph: dict) -> dict:
    return {x["_gvid"]: x for x in graph["objects"] if "pos" in x}


def edge_name(graph_nodes: dict, edge: dict) -> str:
    return f"{graph_nodes[edge['tail']]['name']}→{graph_nodes[edge['head']]['name']}"


def bezier_polyline(points: list, steps: int = 16) -> list:
    """sample a piecewise cubic Bézier into a polyline"""
    out = [tuple(points[0])]
    for i in range(0, len(points) - 3, 3):
        p0, p1, p2, p3 = points[i : i + 4]
        for k in range(1, steps + 1):
            t = k / steps
            u = 1.0 - t
            out.append(
                tuple(
                    u**3 * a + 3 * u * u * t * b + 3 * u * t * t * c + t**3 * d
                    for a, b, c, d in zip(p0, p1, p2, p3)
                )
            )
    return out


def crossings(p: list, q: list) -> list:
    """proper intersection points of two polylines"""

    def side(a, b, c) -> int:
        v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        return (v > 1e-9) - (v < -1e-9)

    out = []
    for a, b in zip(p, p[1:]):
        for c, d in zip(q, q[1:]):
            if side(a, b, c) * side(a, b, d) < 0 and side(c, d, a) * side(c, d, b) < 0:
                out.append(((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0))
    return out


def near(box: tuple[float, float, float, float], point: tuple) -> bool:
    left, bottom, right, top = box
    x, y = point
    return (
        left - _SLACK <= x <= right + _SLACK and bottom - _SLACK <= y <= top + _SLACK
    )


def label_extent(edge: dict) -> tuple[float, float, float, float]:
    """edge label box as (left, bottom, right, top) over all its text lines"""
    size = next(x["size"] for x in edge["_ldraw_"] if x["op"] == "F")
    lines = [x for x in edge["_ldraw_"] if x["op"] == "T"]
    return (
        min(x["pt"][0] - x["width"] / 2.0 for x in lines),
        min(x["pt"][1] for x in lines),
        max(x["pt"][0] + x["width"] / 2.0 for x in lines),
        max(x["pt"][1] for x in lines) + size,
    )


def assert_points_down(graph: dict) -> None:
    graph_nodes = nodes(graph)
    bad = [
        edge_name(graph_nodes, e)
        for e in graph["edges"]
        if e.get("constraint") != "false"
        and center_y(graph_nodes[e["head"]]) >= center_y(graph_nodes[e["tail"]])
    ]
    assert not bad, f"constraint edges not pointing down: {bad}"


@pytest.mark.parametrize(
    "newrank",
    (
        pytest.param(
            False,
            marks=pytest.mark.xfail(
                strict=True,
                reason="newrank=false ranks clusters recursively; inter-cluster "
                "edges are soft (class1.c interclust1)",
            ),
        ),
        True,
    ),
    ids=("default", "newrank"),
)
def test_cluster_chain_points_down(newrank: bool):
    """
    A chain ping-ponging between clusters keeps every constraint edge pointing down.
    """

    assert_points_down(load_layout(COLLAPSE, newrank))


def test_sample_points_down():
    """
    Every constraint edge of the sample points down.
    """

    assert_points_down(load_layout(SAMPLE))


@pytest.mark.parametrize(
    "fixture,newrank",
    ((COLLAPSE, True), (SAMPLE, False)),
    ids=("collapse-newrank", "sample"),
)
def test_adjacent_edges_do_not_cross(fixture: str, newrank: bool):
    """
    Two edges sharing a node do not cross away from that node.
    """

    graph = load_layout(fixture, newrank)
    graph_nodes = nodes(graph)
    edges = [
        (
            edge_name(graph_nodes, e),
            {e["tail"], e["head"]},
            bezier_polyline(draw_points(e)),
        )
        for e in graph["edges"]
    ]

    bad = [
        f"{name1} / {name2}"
        for (name1, ends1, p1), (name2, ends2, p2) in combinations(edges, 2)
        if ends1 & ends2
        and any(
            not any(near(bbox(graph_nodes[n]), x) for n in ends1 & ends2)
            for x in crossings(p1, p2)
        )
    ]
    assert not bad, f"adjacent edges cross: {bad}"


@pytest.mark.parametrize("fixture", ("newrank_fill_order.dot", SAMPLE))
def test_edge_labels_do_not_overlap(fixture: str):
    """
    No two edge labels overlap.
    """

    graph = load_layout(fixture)
    graph_nodes = nodes(graph)
    labels = [
        (edge_name(graph_nodes, e), label_extent(e))
        for e in graph["edges"]
        if "_ldraw_" in e
    ]

    bad = [
        f"{name1} / {name2}"
        for (name1, a), (name2, b) in combinations(labels, 2)
        if a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]
    ]
    assert not bad, f"edge labels overlap: {bad}"
