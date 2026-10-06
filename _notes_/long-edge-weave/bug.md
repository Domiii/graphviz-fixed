# Long constraint=false edge

## HEAD

- These fixtures have flat edges. #2368’s `valid_rank` bail-out skips the `flat_reorder` rewrite, so order is left alone
- Notes often sit on the wrong side of Mid
- Bot→Top is straight in the Mid|Note gap
- This is not the 2.43 far-outer weave

## Far-outer dump (conditional)

- True only when `flat_reorder` rewrites order
- Not the explanation of current HEAD on the edge_shapes fixtures

## Rejected

- HoldAndReinsert reintroduces [Graphviz #2368](https://gitlab.com/graphviz/graphviz/-/issues/2368) (routesplines illegal values)
- Slot-fill + `ED_xpenalty` is blocked by #2368 and is not in the tree

## Rewrite / 2.43 case

- Witness: [`fixtures/clean-flats.dot`](fixtures/clean-flats.dot) → [`bug-shown.svg`](bug-shown.svg), numbers in [`bug-numbers.txt`](bug-numbers.txt)
- Back-edge `SpineBot → SpineTop` (`constraint=false`) L/R/L swings; open lane right of spine unused
- Length / turn / horiz-revs / crossings (red vs green corridor): `1885` vs `942` pt · `887°` vs `180°` · `5` vs `1` · `3` vs `0`
- When `flat_reorder` rewrites, it dumps free long-edge `VIRTUAL` onto the far outer of the flat-Note block (left of Notes)
- Adjacent spine rank keeps that `VIRTUAL` right of Spine → side flip; `ED_xpenalty=0` so mincross will not undo it
- Note-rank `maximal_bbox` is the left sliver; green vertical (right of Spine) lies outside → `Pshortestpath` never sees the corridor
