"""
Labeled fan-out edge detour tests.

Labeled edges from one tail to heads on the next rank (label rank in between, no
real nodes on it) must not swing out to a far label virtual node and back. Here
Src→Right sweeps left around its own label, then right over the top of the red
node Below; Src→Below loops far left before reaching the node right below Src;
Src→Far (several ranks down) swings left past both.
"""

# Do not overfit: assert only the obvious, only-correct behavior (an unobstructed
# fan-out edge neither backtracks horizontally nor leaves its endpoints' band).

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.append(os.path.dirname(__file__))
from gvtest import dot  # pylint: disable=wrong-import-position
from test_edge_shapes import (  # pylint: disable=wrong-import-position
    bbox,
    draw_points,
    edge_by_ends,
    node_by_name,
)

FIXTURES = Path(__file__).parent / "label_detour"

_SLACK = 36.0

SRC = "n14"
RIGHT = "n31"  # Src→Right labeled; head right of Src
BELOW = "n21"  # red node right below Src
FAR = "n53"  # Src→Far labeled; head several ranks down


def load_layout(fixture: str) -> dict:
    path = FIXTURES / fixture
    assert path.exists(), f"unexpectedly missing test case {path}"
    return json.loads(dot("json", path))


def edge_points(graph: dict, tail: str, head: str) -> list:
    return draw_points(
        edge_by_ends(graph, node_by_name(graph, tail), node_by_name(graph, head))
    )


def backtrack(points: list) -> float:
    """horizontal travel beyond net horizontal displacement"""
    travel = sum(abs(x1 - x0) for (x0, _), (x1, _) in zip(points, points[1:]))
    return travel - abs(points[-1][0] - points[0][0])


@pytest.mark.parametrize("head", (RIGHT, BELOW, FAR))
def test_fanout_edge_no_backtrack(head: str):
    """
    Src→head does not swing sideways and back.
    """

    graph = load_layout("fanout_label_detour.dot")
    edge = edge_points(graph, SRC, head)

    back = backtrack(edge)
    assert back <= _SLACK, f"{SRC}→{head} backtracks {back:.1f} > {_SLACK}"


@pytest.mark.parametrize("head", (RIGHT, BELOW, FAR))
def test_fanout_edge_near_endpoints(head: str):
    """
    Src→head stays within its endpoints' horizontal band plus slack.
    """

    graph = load_layout("fanout_label_detour.dot")
    tail_node = node_by_name(graph, SRC)
    head_node = node_by_name(graph, head)
    edge = edge_points(graph, SRC, head)

    low = min(bbox(tail_node)[0], bbox(head_node)[0]) - _SLACK
    high = max(bbox(tail_node)[2], bbox(head_node)[2]) + _SLACK
    xs = [x for x, _ in edge]
    assert low <= min(xs), f"{SRC}→{head} reaches x={min(xs):.1f} < {low:.1f}"
    assert max(xs) <= high, f"{SRC}→{head} reaches x={max(xs):.1f} > {high:.1f}"
