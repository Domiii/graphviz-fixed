"""
Labelled flat edges between clusters: the label node takes a slot the clusters allow.
"""

import os
import subprocess
import sys
from pathlib import Path

sys.path.append(os.path.dirname(__file__))
from gvtest import which  # pylint: disable=wrong-import-position


def test_flat_label_cluster():
    """
    the label node of a flat edge leaving a cluster is not slotted on the cluster's far
    side (x positioning then failed with “trouble in init_rank”)
    """
    input = Path(__file__).parent / "flat_label_cluster.dot"
    assert input.exists(), "unexpectedly missing test case"
    proc = subprocess.run(
        [which("dot"), "-Tplain", input], capture_output=True, text=True, check=False
    )
    assert proc.returncode == 0, proc.stderr
    assert "trouble in init_rank" not in proc.stderr
