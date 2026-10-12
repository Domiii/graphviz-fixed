# Local custom build

* For graphviz-fixed work, **always** use the repo-local install — never system `/tmp` prefixes (`/tmp/gvp`, `mktemp`, etc.; wiped aggressively).
  * Repo: `$GRAPHVIZ_FIXED`
  * Prefix: `$GRAPHVIZ_FIXED/dist`
* Default env before any `dot` for this tree:
  ```
  export GRAPHVIZ_PREFIX=$GRAPHVIZ_FIXED/dist
  export PATH="$GRAPHVIZ_PREFIX/bin:$PATH"
  export LD_LIBRARY_PATH="$GRAPHVIZ_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
  ```
* Stock Graphviz only when the user explicitly wants default: `/usr/bin/dot`.
* Configure / rebuild / install: `dev/build.sh` (`JOBS`, default 4).
* Test: `dev/test.sh <test file or node id> [pytest args]`: targeted, serial, against `dist`; never the full suite.
* Comparisons and regression checks are pytest cases, not ad-hoc scripts.
