#!/usr/bin/env bash
# Rebuild (dev/build.sh), then one targeted, serial pytest selection of
# tests/test_port_fan.py; prints warnings, errors and failures only.
# Usage: dev/retest.sh [-k expr]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
"${ROOT}/dev/build.sh" 2>&1 | rg "warning|error|built:" || true
PYTHON=/usr/bin/python3 "${ROOT}/dev/test.sh" tests/test_port_fan.py -q \
  -p no:cacheprovider --tb=line "$@" 2>&1 |
  rg "^E |^FAILED|passed|failed|big fan" | cut -c1-400 || true
