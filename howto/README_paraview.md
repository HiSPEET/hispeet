## Compare datasets

Let be given datasets `Data1` and `Data2` both containing variable `v` with identical shape. This is how to compare them.

1. Create `Calculator1` by pressing the calculator symbol with `Data1` activated.
2. In the `Properties` menu of `Calculator1` create `v1 = v`.
3. Analogously, create `Calculator2` with `v2 = v` to dataset `Data2`.
4. Activate & mark both calculators and connect them by choosing `Append Attributes`  from menu `Filters/Alphabetical`. This creates `AppendAttributes1`.
5. Create `Calculator3` by pressing the calculator symbol with `AppendAttributes1` activated.
6. In the `Properties` menu of `AppendAttributes1` create the new variable `dv = v2 - v1`
7. Explore this new variable as usual.