1.  On each level mark active elements
    -  `1000` for refinement,    set with `element % MarkForRefinement()`  
    -  `   0` for retention,     set with `element % MarkForRetention()` 
    -  `  -1` for removal,       set with `element % MarkForRemoval()`

2. Restriction of child marks
    -  set cluster mark `cm = m + o` where `m` is the maximum of corresponding
       child marks and `c` is bitwise true for every child marked for refinement
    -  transfer child cluster marks to parent 
    -  if `cm ≥ 0` upgrade parent marks to `pm = max(0,pm) + 1000 + cm`

3.  Globalize adaptation pattern
    -  copy adaptation mark to ghosts
    -  upgrade element marks according to neighbor refinement
       *  `max(mark, 0)` if neighbor neighbor mark `nm ≥ 1000`
       *  `max(mark, 1000)` if abutting to any child marked for refinement
    -  upgrade marks to the nearest thousand, i.e. add level if any child is refined
