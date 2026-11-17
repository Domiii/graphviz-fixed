"""
`-D(GNEA)name=val`: command-line attribute defaults that the input graph overrides.
"""

import json
import os
import sys

sys.path.append(os.path.dirname(__file__))
from gvtest import run, which  # pylint: disable=wrong-import-position


def test_soft_defaults():
    """
    unset attributes get the `-D…` value; any value set in the input wins
    """
    source = """digraph {
      label=file
      node [fontcolor=red]
      edge [color=blue]
      a [color=green]
      a -> b
      subgraph s { node [fontcolor=yellow]; c }
    }"""
    flags = [
        "-DGlabel=soft",
        "-DGbgcolor=black",
        "-DNfontcolor=white",
        "-DNcolor=gray",
        "-DEcolor=gray",
        "-DEfontcolor=white",
    ]
    out = json.loads(run(which("dot"), "-Tjson0", *flags, input=source))
    nodes = {o["name"]: o for o in out["objects"] if "nodes" not in o}
    (edge,) = out["edges"]

    # input wins: graph attribute, root/subgraph defaults, object attribute
    assert out["label"] == "file"
    assert nodes["a"]["fontcolor"] == "red"
    assert nodes["c"]["fontcolor"] == "yellow"
    assert nodes["a"]["color"] == "green"
    assert edge["color"] == "blue"

    # unset in the input: soft default applies
    assert out["bgcolor"] == "black"
    assert nodes["b"]["color"] == "gray"
    assert nodes["c"]["color"] == "gray"
    assert edge["fontcolor"] == "white"
