#!/usr/bin/env bash
# Configure (once), build and install this tree into $ROOT/dist.
# JOBS: build parallelism; default 4.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
JOBS="${JOBS:-4}"

cd "${ROOT}"
if [ ! -f build/CMakeCache.txt ]; then
  nice -n 10 cmake -DCMAKE_INSTALL_PREFIX="${ROOT}/dist" -B build -S .
fi
nice -n 10 cmake --build build -j "${JOBS}"
cmake --install build >/dev/null
LD_LIBRARY_PATH="${ROOT}/dist/lib" "${ROOT}/dist/bin/dot" -c
echo "built: ${ROOT}/dist/bin/dot"
