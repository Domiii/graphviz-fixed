"""
`rankpos`: nodes are pulled toward their `rankpos` differences along the rank.
"""

import os
import sys

import pytest

sys.path.append(os.path.dirname(__file__))
from gvtest import run, which  # pylint: disable=wrong-import-position


def rank_positions(source: str, rankdir: str) -> dict[str, float]:
    """
    node positions along their rank (points): left to right, or top to bottom
    for LR/RL
    """
    out = run(
        which("dot"), "-Tplain", input=source.replace("{", f"{{rankdir={rankdir};", 1)
    )
    pos = {}
    for line in out.splitlines():
        fields = line.split()
        if fields[0] == "node":
            pos[fields[1]] = 72 * (
                float(fields[2]) if rankdir in ("TB", "BT") else -float(fields[3])
            )
    return pos


@pytest.mark.parametrize("rankdir", ["TB", "BT", "LR", "RL"])
def test_rankpos(rankdir: str):
    """
    hinted nodes take their `rankpos` differences against edge straightening;
    the drawing may shift
    """
    source = (
        "digraph { rankposweight=2; r -> {a b c}; a [rankpos=1000]; b [rankpos=1400] }"
    )
    pos = rank_positions(source, rankdir)
    assert pos["b"] - pos["a"] == pytest.approx(400, abs=1)
    free = rank_positions(source.replace("rankposweight=2", "rankposweight=0"), rankdir)
    assert free["b"] - free["a"] < 200


def test_rankpos_unset():
    """an empty `rankpos` is no hint"""
    source = "digraph { r -> {a b c}; a -> c; b -> a }"
    declared = source.replace("{", '{ node [rankpos=""];', 1)
    assert run(which("dot"), "-Tplain", input=declared) == run(
        which("dot"), "-Tplain", input=source
    )
