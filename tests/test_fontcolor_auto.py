"""
`fontcolor=auto`: black or white text, whichever contrasts more with the label’s background.
"""

import json
import os
import sys

sys.path.append(os.path.dirname(__file__))
from gvtest import run, which  # pylint: disable=wrong-import-position


def test_fontcolor_auto():
    """
    node fill, else enclosing filled cluster, else `pagecolor`; explicit colors are kept
    """
    source = """digraph {
      pagecolor="#111111"
      node [fontcolor=auto]
      light [style=filled fillcolor=yellow]
      dark [style=filled fillcolor=navy]
      plain
      explicit [style=filled fillcolor=navy fontcolor=red]
      subgraph cluster_c {
        label=c fontcolor=auto style=filled fillcolor=lightyellow
        inner
      }
    }"""
    out = json.loads(run(which("dot"), "-Tjson", input=source))

    def text_color(obj: dict) -> str:
        (color,) = [op["color"] for op in obj["_ldraw_"] if op["op"] == "c"]
        return color

    colors = {o["name"]: text_color(o) for o in out["objects"]}

    assert colors["light"] == "#000000"
    assert colors["dark"] == "#ffffff"
    assert colors["plain"] == "#ffffff"
    assert colors["inner"] == "#000000"
    assert colors["cluster_c"] == "#000000"
    assert colors["explicit"] == "#ff0000"
