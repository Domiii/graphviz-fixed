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
* Configure / rebuild / install:
  ```
  cd $GRAPHVIZ_FIXED
  cmake -DCMAKE_INSTALL_PREFIX=$PWD/dist -B build -S .
  cmake --build build -j16
  cmake --install build
  LD_LIBRARY_PATH=$GRAPHVIZ_PREFIX/lib $GRAPHVIZ_PREFIX/bin/dot -c
  ```
