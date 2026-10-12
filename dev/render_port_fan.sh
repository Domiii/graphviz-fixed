#!/usr/bin/env bash
# Render pr-artifacts/port-fan/<fixture>.{before,after}.svg from tests/port_fan/*.dot.
# before: BEFORE_PREFIX (a dist built at the base commit); after: this tree's dist.
# Usage: BEFORE_PREFIX=/path/to/base/dist dev/render_port_fan.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${ROOT}/pr-artifacts/port-fan"
: "${BEFORE_PREFIX:?set BEFORE_PREFIX to a dist built at the base commit}"
mkdir -p "${OUT}"
for f in "${ROOT}"/tests/port_fan/*.dot; do
  name="$(basename "${f}" .dot)"
  name="${name//_/-}"
  LD_LIBRARY_PATH="${BEFORE_PREFIX}/lib" nice -n 10 "${BEFORE_PREFIX}/bin/dot" -Tsvg "${f}" -o "${OUT}/${name}.before.svg" 2>/dev/null || true
  LD_LIBRARY_PATH="${ROOT}/dist/lib" nice -n 10 "${ROOT}/dist/bin/dot" -Tsvg "${f}" -o "${OUT}/${name}.after.svg"
  echo "wrote ${name}"
done
