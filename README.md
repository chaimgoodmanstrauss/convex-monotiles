# Convex Monotiles

Reinhardt-style families of convex monotiles: triangle, quadrilateral, 15 pentagon types, 3 hexagon types.
Labels: family · k (isohedral number: tile orbits under the tiling's symmetry) · orbifold name of the symmetry group.

- `families.js` — the defining equations of each family (Wikipedia "Pentagonal tiling" conditions; vertices CCW, side a = E→A) and vertex-drag projection.
- `tiler.js` — periodic tilings as (wallpaper group, k tile orbits); overlap-area search, contact equations for continuation while dragging.
- `solutions.js` — generated offline: a member of each family plus a perfect tiling of it.
- `edgesym-derive.mjs` → `edgesym-data.js` (experiment, run `node edgesym-derive.mjs`): from each stored tiling, the linear constraints its edge symmetries impose — corner-angle sums at each vertex point, edge-length sums along each line — plus a gluing tree that rebuilds the cell from the tile.
- `edgesym.js` — parametrizes each type directly by those constraints (angles affine, sides linear, closure projected linearly; one scalar closure equation when only two side lengths are free), lays the tile down by turtle, and glues the cell together; no projection onto family equations and no tiling solver. Toggle in the app: "constraints: edge symmetries / family equations".
