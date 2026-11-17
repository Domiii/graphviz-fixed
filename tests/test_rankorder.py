"""
`rankorder`: nodes (and clusters, via their members) keep ascending `rankorder` within a rank.
"""

import os
import sys

import pytest

sys.path.append(os.path.dirname(__file__))
from gvtest import run, which  # pylint: disable=wrong-import-position


def drawing_order(source: str, rankdir: str, names: set[str]) -> list[str]:
    """`names` in the order their rank is drawn: left to right, or top to bottom for LR/RL"""
    out = run(which("dot"), "-Tplain", input=source.replace("{", f"{{rankdir={rankdir};", 1))
    pos = {}
    for line in out.splitlines():
        fields = line.split()
        if fields[0] == "node" and fields[1] in names:
            pos[fields[1]] = float(fields[2]) if rankdir in ("TB", "BT") else -float(fields[3])
    return sorted(pos, key=pos.get)


@pytest.mark.parametrize("rankdir", ["TB", "BT", "LR", "RL"])
def test_rankorder(rankdir: str):
    """
    hinted nodes and clusters keep their hint order against crossings and input order;
    unhinted nodes are placed freely
    """
    crossing = "digraph { r -> {a b}; a -> x; b -> y; x [rankorder=2]; y [rankorder=1] }"
    assert drawing_order(crossing, rankdir, {"x", "y"}) == ["y", "x"]
    assert drawing_order(crossing, rankdir, {"a", "b"}) == ["b", "a"]

    clusters = """digraph { r -> {p q}
      subgraph cluster_A { p -> {p1 p2} } subgraph cluster_B { q -> {q1 q2} }
      p1 [rankorder=5] q1 [rankorder=1] }"""
    assert drawing_order(clusters, rankdir, {"p", "q"}) == ["q", "p"]


def test_rankorder_unset():
    """an empty `rankorder` is no hint"""
    source = "digraph { r -> {a b c}; a -> c; b -> a }"
    declared = source.replace("{", '{ node [rankorder=""];', 1)
    assert run(which("dot"), "-Tplain", input=declared) == run(
        which("dot"), "-Tplain", input=source
    )
