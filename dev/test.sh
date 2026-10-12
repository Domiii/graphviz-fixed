#!/usr/bin/env bash
# Targeted pytest against this tree's $ROOT/dist. Serial, niced.
# Usage: dev/test.sh tests/test_port_fan.py [-k expr] [pytest args...]
# Refuses to run without a target: never the full suite.
# PYTHON: interpreter that has pytest; default python3.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON="${PYTHON:-python3}"

if [ "$#" -eq 0 ]; then
  echo "usage: $0 <test file or node id> [pytest args...]" >&2
  exit 2
fi
if [ ! -x "${ROOT}/dist/bin/dot" ]; then
  echo "no ${ROOT}/dist/bin/dot; run dev/build.sh" >&2
  exit 1
fi

cd "${ROOT}"
export PATH="${ROOT}/dist/bin:${PATH}"
export LD_LIBRARY_PATH="${ROOT}/dist/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
exec nice -n 10 "${PYTHON}" -m pytest -p no:xdist "$@"
