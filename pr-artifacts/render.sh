#!/usr/bin/env bash
# Render before/after SVGs for tests/edge_shapes/*.dot
# before = system /usr/bin/dot (unfixed)
# after  = local fixed install (GRAPHVIZ_PREFIX or /tmp/graphviz-prefix.*)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$(cd "$(dirname "$0")" && pwd)/edge_shapes"
FIXTURES="${ROOT}/tests/edge_shapes"

BEFORE_DOT="${BEFORE_DOT:-/usr/bin/dot}"

if [ -z "${GRAPHVIZ_PREFIX:-}" ]; then
  GRAPHVIZ_PREFIX="$(ls -d /tmp/graphviz-prefix.* 2>/dev/null | head -1 || true)"
fi

if [ -z "${GRAPHVIZ_PREFIX}" ] || [ ! -x "${GRAPHVIZ_PREFIX}/bin/dot" ]; then
  echo "No fixed install found. Build per DEVELOPERS.md, e.g.:" >&2
  echo "  PREFIX=\$(mktemp -d -t graphviz-prefix.XXXXXX)" >&2
  echo "  cmake -DCMAKE_INSTALL_PREFIX=\${PREFIX} -B build -S ." >&2
  echo "  cmake --build build && cmake --install build" >&2
  echo "  GRAPHVIZ_PREFIX=\${PREFIX} $0" >&2
  exit 1
fi

AFTER_DOT="${GRAPHVIZ_PREFIX}/bin/dot"

mkdir -p "${OUT}"

echo "BEFORE_DOT=${BEFORE_DOT}"
"${BEFORE_DOT}" -V
echo "AFTER_DOT=${AFTER_DOT}"
env PATH="${GRAPHVIZ_PREFIX}/bin:${PATH}" \
  LD_LIBRARY_PATH="${GRAPHVIZ_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}" \
  "${AFTER_DOT}" -V

for dotfile in "${FIXTURES}"/*.dot; do
  name="$(basename "${dotfile}" .dot)"
  "${BEFORE_DOT}" -Tsvg "${dotfile}" -o "${OUT}/${name}.before.svg"
  env PATH="${GRAPHVIZ_PREFIX}/bin:${PATH}" \
    LD_LIBRARY_PATH="${GRAPHVIZ_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}" \
    "${AFTER_DOT}" -Tsvg "${dotfile}" -o "${OUT}/${name}.after.svg"
  echo "wrote ${name}.before.svg ${name}.after.svg"
done
