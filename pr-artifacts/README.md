before = system `/usr/bin/dot` (unfixed). after = local fixed build (`GRAPHVIZ_PREFIX` or `dist`).

- `edge_shapes/`: re-run `GRAPHVIZ_PREFIX=$GRAPHVIZ_FIXED/dist ./pr-artifacts/render.sh`
- `complex-weave/`: render from `tests/complex_weave/*.dot` and `_notes_/complex-weave/sample2.dot` (see PR notes)
- `label-detour/`: render from `tests/label_detour/*.dot`; before = local `dist` build at HEAD (bug absent in system 2.43)
- `many-parallel/`: render from `tests/many_parallel/*.dot`; before = local `dist` build at HEAD (`tasks_and_containers_coarse.before` = sample without `newrank`, `.newrank-only` = with `newrank` at HEAD; `cluster_rank_collapse.after` adds `newrank=true`)
