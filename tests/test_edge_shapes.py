"""
Long constraint=false edge far-outer layout tests.

A free long back-edge must sit on Mid's free side opposite the Note — not in
the Mid|Note gap, and not past the Note's far outer — when that opposite side
is free. Flip fixtures only forbid wrapping both Notes' far outers.

Corridor and Note/Mid LR pins use intended geometry from the .dot (flat edge
direction / rank=same listing), not landed placement.
"""

# Do not overfit: assert only the obvious, only-correct behavior (e.g. never on a
# Note's far outer when the opposite side is free and crossing-free). Where several
# correct corridors exist, assert invariants that rule out the weave, not one polyline.

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.append(os.path.dirname(__file__))
from gvtest import dot  # pylint: disable=wrong-import-position

FIXTURES = Path(__file__).parent / "edge_shapes"

_POINT_OPS = frozenset({"p", "P", "b", "B"})
_POINTS_PER_INCH = 72.0


def node_by_name(graph: dict, name: str) -> dict:
    matches = [x for x in graph["objects"] if x.get("name") == name]
    assert len(matches) == 1, f"expected one node named {name!r}"
    return matches[0]


def edge_by_ends(graph: dict, tail: dict, head: dict) -> dict:
    matches = [
        x
        for x in graph["edges"]
        if x["tail"] == tail["_gvid"] and x["head"] == head["_gvid"]
    ]
    assert len(matches) == 1, "expected one edge between the given endpoints"
    return matches[0]


def draw_points(obj: dict) -> list:
    for item in obj.get("_draw_", []):
        if item.get("op") in _POINT_OPS and "points" in item:
            return item["points"]
    raise AssertionError("no polygon/polyline points in _draw_")


def bbox(node: dict) -> tuple[float, float, float, float]:
    """node box as (left, bottom, right, top) from JSON pos/width/height"""
    cx, cy = (float(p) for p in node["pos"].split(","))
    half_w = float(node["width"]) * _POINTS_PER_INCH / 2.0
    half_h = float(node["height"]) * _POINTS_PER_INCH / 2.0
    return cx - half_w, cy - half_h, cx + half_w, cy + half_h


def center_y(node: dict) -> float:
    return float(node["pos"].split(",")[1])


def assert_left_of(left: dict, right: dict) -> None:
    assert (
        bbox(left)[0] < bbox(right)[0]
    ), f"expected {left.get('name')} left of {right.get('name')}"


def x_at_y(points: list, y: float) -> list[float]:
    """interpolate polyline x where segments cross horizontal line y"""
    xs = []
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 == y1:
            if y0 == y:
                xs.extend([x0, x1])
            continue
        if (y0 - y) * (y1 - y) > 0:
            continue
        t = (y - y0) / (y1 - y0)
        xs.append(x0 + t * (x1 - x0))
    assert xs, f"polyline does not cross y={y}"
    return xs


def load_layout(fixture: str) -> dict:
    path = FIXTURES / fixture
    assert path.exists(), f"unexpectedly missing test case {path}"
    return json.loads(dot("json", path))


def back_edge_points(graph: dict) -> list:
    bot = node_by_name(graph, "Bot")
    top = node_by_name(graph, "Top")
    return draw_points(edge_by_ends(graph, bot, top))


def far_outer(note: dict, mid: dict) -> tuple[str, float]:
    """side of the note away from mid, as ('L' or 'R', x)"""
    note_left, _, note_right, _ = bbox(note)
    mid_left, _, mid_right, _ = bbox(mid)
    if (note_left + note_right) < (mid_left + mid_right):
        return "L", note_left
    return "R", note_right


def all_beyond_far_outer(
    edge_points: list, y_node: dict, note: dict, mid: dict
) -> bool:
    side, limit = far_outer(note, mid)
    xs = x_at_y(edge_points, center_y(y_node))
    if side == "L":
        return all(x < limit for x in xs)
    return all(x > limit for x in xs)


def assert_not_beyond_far_outer(edge_points: list, mid: dict, note: dict) -> None:
    """Bot→Top at mid's center y must not be strictly past the note's far outer"""
    side, limit = far_outer(note, mid)
    xs = x_at_y(edge_points, center_y(mid))
    name = note.get("name", "?")
    if side == "L":
        ok = all(x >= limit for x in xs)
    else:
        ok = all(x <= limit for x in xs)
    assert ok, f"Bot→Top beyond {name} far outer"


def assert_on_free_side(edge_points: list, mid: dict, free_side: str) -> None:
    """
    At Mid's center_y, Bot→Top must be on Mid's free side ('L' or 'R').
    free_side is intended geometry, not inferred from landed Note placement.
    """
    mid_left, _, mid_right, _ = bbox(mid)
    xs = x_at_y(edge_points, center_y(mid))
    if free_side == "R":
        ok = all(x >= mid_right for x in xs)
        where = "right of Mid"
    else:
        assert free_side == "L", f"free_side must be 'L' or 'R', got {free_side!r}"
        ok = all(x <= mid_left for x in xs)
        where = "left of Mid"
    assert ok, f"Bot→Top not on free side ({where})"


@pytest.mark.parametrize(
    "fixture,note_left_of_mid,free_side",
    (
        ("right_stays_right.dot", True, "R"),
        ("left_stays_left.dot", False, "L"),
    ),
)
def test_not_beyond_note_far_outer(
    fixture: str, note_left_of_mid: bool, free_side: str
):
    """
    Single Note beside Mid: pin intended Note|Mid LR; Bot→Top on free side.
    """

    graph = load_layout(fixture)
    mid = node_by_name(graph, "Mid")
    note = node_by_name(graph, "Note")
    edge = back_edge_points(graph)

    if note_left_of_mid:
        assert_left_of(note, mid)
    else:
        assert_left_of(mid, note)

    assert_not_beyond_far_outer(edge, mid, note)
    assert_on_free_side(edge, mid, free_side)


@pytest.mark.parametrize(
    "fixture",
    (
        "flip_free_R_then_L.dot",
        "flip_free_L_then_R.dot",
    ),
)
def test_flip_not_full_width_weave(fixture: str):
    """
    Notes at MidA and MidB: pin intended Note sides; forbid wrapping both far outers.
    """

    graph = load_layout(fixture)
    mid_a = node_by_name(graph, "MidA")
    mid_b = node_by_name(graph, "MidB")
    note_l = node_by_name(graph, "NoteL")
    note_r = node_by_name(graph, "NoteR")
    edge = back_edge_points(graph)

    if fixture == "flip_free_R_then_L.dot":
        assert_left_of(note_l, mid_a)
        assert_left_of(mid_b, note_r)
    else:
        assert_left_of(mid_a, note_r)
        assert_left_of(note_l, mid_b)

    def rank_mid(note: dict) -> dict:
        if abs(center_y(note) - center_y(mid_a)) <= abs(
            center_y(note) - center_y(mid_b)
        ):
            return mid_a
        return mid_b

    outside_l = all_beyond_far_outer(edge, note_l, note_l, rank_mid(note_l))
    outside_r = all_beyond_far_outer(edge, note_r, note_r, rank_mid(note_r))
    assert not (
        outside_l and outside_r
    ), "Bot→Top wraps both Notes' far outers (full-width weave)"


def test_not_notes_far_outer():
    """
    LeftNote left of Mid, RightBlock right: Bot→Top on Mid's free side (right).
    """

    graph = load_layout("not_notes_far_outer.dot")
    mid = node_by_name(graph, "Mid")
    right = node_by_name(graph, "RightBlock")
    note = node_by_name(graph, "LeftNote")
    edge = back_edge_points(graph)

    assert_left_of(note, mid)
    assert_left_of(mid, right)
    assert_not_beyond_far_outer(edge, mid, note)
    assert_on_free_side(edge, mid, "R")
