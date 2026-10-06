"""
Labeled constraint=false back-edge weave tests.

A back edge from a low rank to a high rank must take a near corridor, not detour
around the far side of the whole graph. Back edges must not dominate graph width,
and sibling "judged" back edges must stay near their endpoints and off unrelated
edge labels.
"""

# Do not overfit: assert only the obvious, only-correct behavior (e.g. never wrap
# the far outer of a column the edge has no business with). Where several correct
# corridors exist, assert invariants that rule out the weave, not one polyline.

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.append(os.path.dirname(__file__))
from gvtest import dot  # pylint: disable=wrong-import-position
from test_edge_shapes import (  # pylint: disable=wrong-import-position
    bbox,
    center_y,
    draw_points,
    edge_by_ends,
    node_by_name,
    x_at_y,
)

FIXTURES = Path(__file__).parent / "complex_weave"

_SLACK = 36.0
_EXCURSION_SLACK = 216.0


def load_layout(fixture: str) -> dict:
    path = FIXTURES / fixture
    assert path.exists(), f"unexpectedly missing test case {path}"
    return json.loads(dot("json", path))


def nodes(graph: dict) -> list:
    return [x for x in graph["objects"] if "pos" in x]


def edge_points(graph: dict, tail: str, head: str) -> list:
    return draw_points(
        edge_by_ends(graph, node_by_name(graph, tail), node_by_name(graph, head))
    )


def label_box(edge: dict) -> tuple[float, float, float, float]:
    """edge label box as (left, bottom, right, top) from _ldraw_ font and text ops"""
    size = next(x["size"] for x in edge["_ldraw_"] if x["op"] == "F")
    text = next(x for x in edge["_ldraw_"] if x["op"] == "T")
    assert text["align"] == "c", "expected a centered edge label"
    x, y = text["pt"]
    half_w = text["width"] / 2.0
    return x - half_w, y - size / 2.0, x + half_w, y + size


def band_right(graph: dict, tail: dict, head: dict) -> float:
    """rightmost node right edge among nodes centered within the tail..head y-span"""
    lo, hi = sorted((center_y(tail), center_y(head)))
    return max(bbox(n)[2] for n in nodes(graph) if lo <= center_y(n) <= hi)


def content_right(graph: dict) -> float:
    """rightmost node, constraint-edge spline, or constraint-edge label"""
    rights = [bbox(n)[2] for n in nodes(graph)]
    for edge in graph["edges"]:
        if edge.get("constraint") == "false":
            continue
        rights.extend(x for x, _ in draw_points(edge))
        if "_ldraw_" in edge:
            rights.append(label_box(edge)[2])
    return max(rights)


def test_back_edge_not_around_far_column():
    """
    Low (bottom of one column) back to Top must not wrap the other column's far outer.
    """

    graph = load_layout("far_side_back_edge.dot")
    low = node_by_name(graph, "Low")
    edge = edge_points(graph, "Low", "Top")
    column = [node_by_name(graph, n) for n in ("R1", "R2", "R3")]

    far_right = bbox(column[0])[0] > bbox(low)[0]
    for node in column:
        xs = x_at_y(edge, center_y(node))
        left, _, right, _ = bbox(node)
        if far_right:
            ok = all(x <= right for x in xs)
        else:
            ok = all(x >= left for x in xs)
        assert ok, f"Low→Top wraps {node['name']} far outer"


@pytest.mark.parametrize(
    "tail",
    ("evaluateChecks", "runRequirementsJudge"),
)
def test_back_edge_excursion_bound(tail: str):
    """
    Back edge to ToolRunner.run stays within the nodes of its y-span plus slack.
    """

    graph = load_layout("sample2_min.dot")
    tail_node = node_by_name(graph, tail)
    head_node = node_by_name(graph, "ToolRunner.run")
    edge = edge_points(graph, tail, "ToolRunner.run")

    limit = band_right(graph, tail_node, head_node) + _SLACK
    reach = max(x for x, _ in edge)
    assert reach <= limit, f"{tail}→ToolRunner.run reaches x={reach} > {limit}"


def test_width_not_dominated_by_back_edges():
    """
    Graph bbox extends at most slack past nodes and constraint edges with labels.
    """

    graph = load_layout("sample2_min.dot")
    right = float(graph["bb"].split(",")[2])
    limit = content_right(graph) + _SLACK
    assert right <= limit, f"graph right x={right} > {limit}"


@pytest.mark.parametrize(
    "tail",
    ("resolve_email", "resolve_missing_email"),
)
def test_judged_edge_near_endpoints(tail: str):
    """
    Judged back edge to runEmailStep stays in its endpoints' horizontal band.
    """

    graph = load_layout("sample2_min.dot")
    tail_node = node_by_name(graph, tail)
    head_node = node_by_name(graph, "runEmailStep")
    edge = edge_points(graph, tail, "runEmailStep")

    low = min(bbox(tail_node)[0], bbox(head_node)[0]) - _EXCURSION_SLACK
    high = max(bbox(tail_node)[2], bbox(head_node)[2]) + _EXCURSION_SLACK
    xs = [x for x, _ in edge]
    assert low <= min(xs), f"{tail}→runEmailStep reaches x={min(xs)} < {low}"
    assert max(xs) <= high, f"{tail}→runEmailStep reaches x={max(xs)} > {high}"


def test_judged_edge_off_finish_label():
    """
    resolve_missing_email→runEmailStep must not cross the finishStep→end label.
    """

    graph = load_layout("sample2_min.dot")
    finish = edge_by_ends(
        graph,
        node_by_name(graph, "finishStep"),
        node_by_name(graph, "JourneyRunState.end"),
    )
    left, bottom, right, top = label_box(finish)
    edge = edge_points(graph, "resolve_missing_email", "runEmailStep")

    for y in (bottom, (bottom + top) / 2.0, top):
        xs = x_at_y(edge, y)
        assert not any(
            left <= x <= right for x in xs
        ), f"resolve_missing_email→runEmailStep crosses finish label at y={y}"
