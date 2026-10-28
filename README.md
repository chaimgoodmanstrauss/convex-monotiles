# Convex Monotiles

Reinhardt-style families of convex monotiles: triangle, quadrilateral, 15 pentagon types, 3 hexagon types.
Labels: family · k (isohedral number: tile orbits under the tiling's symmetry) · orbifold name of the symmetry group.

- `families.js` — the defining equations of each family (Wikipedia "Pentagonal tiling" conditions; vertices CCW, side a = E→A) and vertex-drag projection.
- `tiler.js` — periodic tilings as (wallpaper group, k tile orbits); overlap-area search, contact equations for continuation while dragging.
- `solutions.js` — generated offline: a member of each family plus a perfect tiling of it.
