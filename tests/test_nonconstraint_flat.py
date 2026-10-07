"""
constraint=false flat edge tests.

A constraint=false edge between two nodes of one rank imposes no left-to-right
order on them. Here only it could force y left of x and cross a→x with b→y.
"""

# Do not overfit: assert only the obvious, only-correct behavior (the crossing-free
# order of the children follows their parents).

import json
import os
import sys
from pathlib import Path

sys.path.append(os.path.dirname(__file__))
from gvtest import dot  # pylint: disable=wrong-import-position
from test_edge_shapes import (  # pylint: disable=wrong-import-position
    assert_left_of,
    node_by_name,
)

FIXTURES = Path(__file__).parent / "nonconstraint_flat"


def test_flat_imposes_no_order():
    """
    Children keep their parents' order, so a→x and b→y do not cross.
    """

    path = FIXTURES / "flat_no_order.dot"
    assert path.exists(), f"unexpectedly missing test case {path}"
    graph = json.loads(dot("json", path))

    assert_left_of(node_by_name(graph, "a"), node_by_name(graph, "b"))
    assert_left_of(node_by_name(graph, "x"), node_by_name(graph, "y"))
