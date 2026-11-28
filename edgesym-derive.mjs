// Offline: derive, for each family's stored tiling, the linear constraints that its edge symmetries impose
// (corner-angle sums at every vertex point, edge-length sums along every line of edges) and a gluing tree
// that rebuilds the whole cell from the tile alone. Writes edgesym-data.js.   Run: node edgesym-derive.mjs
import {writeFileSync} from 'node:fs';
import {FAMILIES, polyOf} from './families.js';
import {SOLS} from './solutions.js';
import {buildCell, findContacts, solveContacts, GROUPS, polyArea} from './tiler.js';

const TOL = 1e-6, R = 3;
const out = {};

for (const fam of FAMILIES) {
  const sol = SOLS[fam.key]; if (!sol || fam.n < 5) continue;
  const n = fam.n, cfg = {group: fam.group, k: fam.k, flips: sol.flips}, P = polyOf(sol.P);
  if (polyArea(P) <= 0) throw new Error(fam.key + ': tile not CCW');
  let x = sol.x; x = solveContacts(cfg, P, x, findContacts(cfg, P, x, 1e-5)).x;
  const {lat, tiles} = buildCell(cfg, P, x), ord = GROUPS[fam.group].order;
  const detOf = c => { const g = c % ord, j = Math.floor(c / ord); const dg = (fam.group === 'pgg' && g >= 2) ? -1 : 1; return dg * (sol.flips[j] ? -1 : 1); };
  // patch of labeled tiles: W[label] = world position of that corner
  const patch = [];
  for (let I = -R; I <= R; I++) for (let J = -R; J <= R; J++) tiles.forEach((t, c) => {
    const s = detOf(c), sx = I * lat.t1.x + J * lat.t2.x, sy = I * lat.t1.y + J * lat.t2.y;
    const W = Array.from({length: n}, (_, L) => { const v = t.poly[s > 0 ? L : n - 1 - L]; return {x: v.x + sx, y: v.y + sy}; });
    patch.push({c, I, J, det: s, W});
  });
  const idx = (c, I, J) => patch.findIndex(p => p.c === c && p.I === I && p.J === J);
  const central = patch.filter(p => p.I === 0 && p.J === 0);
  const near = (a, b) => Math.hypot(a.x - b.x, a.y - b.y) < TOL;
  const onSegInterior = (p, a, b) => {
    const ex = b.x - a.x, ey = b.y - a.y, l2 = ex * ex + ey * ey, t = ((p.x - a.x) * ex + (p.y - a.y) * ey) / l2;
    return t > 1e-6 && t < 1 - 1e-6 && Math.abs(ex * (p.y - a.y) - ey * (p.x - a.x)) / Math.sqrt(l2) < TOL;
  };

  // ---- angle constraints (units of pi): corners meeting at a point sum to 2, or 1 if the point is inside another tile's edge
  const angRows = [], angWhy = [], seenPts = [];
  for (const T of central) for (let L = 0; L < n; L++) {
    const p = T.W[L]; if (seenPts.some(q => near(p, q))) continue; seenPts.push(p);
    const row = new Array(n).fill(0); let flat = false;
    for (const U of patch) for (let M = 0; M < n; M++) {
      if (near(U.W[M], p)) row[M]++;
      else if (onSegInterior(p, U.W[M], U.W[(M + 1) % n])) flat = true;
    }
    angRows.push([...row, flat ? 1 : 2]);
  }
  angRows.push([...new Array(n).fill(1), n - 2]); // angle sum of the polygon

  // ---- side constraints: along each line of edges, the two sides' edge lengths must reach the same points
  // unknowns: positions of points on the line, then side lengths s_0..s_{n-1}; each edge gives x_end - x_start - (+/-)s = 0
  const sideRows = [];
  const doneLines = [];
  for (const T of central) for (let L = 0; L < n; L++) {
    const a = T.W[L], b = T.W[(L + 1) % n], ux = (b.x - a.x), uy = (b.y - a.y), ul = Math.hypot(ux, uy), dx = ux / ul, dy = uy / ul;
    const onLine = p => Math.abs((p.x - a.x) * dy - (p.y - a.y) * dx) < TOL;
    const key = doneLines.find(o => Math.abs(o.dx * dy - o.dy * dx) < 1e-9 && onLine(o.a));
    if (key) continue; doneLines.push({a, dx, dy});
    const pts = [], edges = [];
    const ptIdx = p => { const c = (p.x - a.x) * dx + (p.y - a.y) * dy; let i = pts.findIndex(q => Math.abs(q - c) < TOL); if (i < 0) { pts.push(c); i = pts.length - 1; } return i; };
    for (const U of patch) for (let M = 0; M < n; M++) {
      const p = U.W[M], q = U.W[(M + 1) % n];
      if (!onLine(p) || !onLine(q)) continue;
      const i = ptIdx(p), j = ptIdx(q);
      edges.push({i, j, M, sgn: pts[j] > pts[i] ? 1 : -1});
    }
    const np = pts.length, rows = edges.map(e => { const r = new Array(np + n).fill(0); r[e.j] += 1; r[e.i] -= 1; r[np + e.M] -= e.sgn; return r; });
    // eliminate the point positions; rows left with zero point-part are relations among the sides
    const red = rref(rows, np);
    for (const r of red) if (r.slice(0, np).every(v => Math.abs(v) < 1e-9) && r.slice(np).some(v => Math.abs(v) > 1e-9)) sideRows.push(r.slice(np));
  }

  // ---- gluing tree: tile B is placed from a placed tile A by a shared corner and a shared edge ray
  const targets = new Set([0, idx(0, 1, 0), idx(0, 0, 1)]);
  const repOf = new Map(); // cell index -> patch index
  const start = idx(0, 0, 0), from = new Map([[start, null]]), queue = [start];
  while (queue.length) {
    const ai = queue.shift(), A = patch[ai];
    if (!repOf.has(A.c)) repOf.set(A.c, ai);
    for (let bi = 0; bi < patch.length; bi++) {
      if (from.has(bi)) continue;
      const B = patch[bi]; let g = null;
      for (let a = 0; a < n && !g; a++) for (let b = 0; b < n && !g; b++) {
        if (!near(A.W[a], B.W[b])) continue;
        for (const da of [1, n - 1]) for (const db of [1, n - 1]) {
          if (g) break;
          const a2 = (a + da) % n, b2 = (b + db) % n;
          const va = {x: A.W[a2].x - A.W[a].x, y: A.W[a2].y - A.W[a].y}, vb = {x: B.W[b2].x - B.W[b].x, y: B.W[b2].y - B.W[b].y};
          const la = Math.hypot(va.x, va.y), lb = Math.hypot(vb.x, vb.y);
          if (Math.abs(va.x * vb.y - va.y * vb.x) / (la * lb) < 1e-7 && va.x * vb.x + va.y * vb.y > 0) g = {from: ai, a, a2, b, b2};
        }
      }
      if (g) { from.set(bi, g); queue.push(bi); }
    }
  }
  const need = new Set([idx(0, 1, 0), idx(0, 0, 1), ...repOf.values()]);
  const steps = [], have = new Set([start]);
  const addPath = t => { if (have.has(t)) return; const g = from.get(t); if (!g) throw new Error(fam.key + ': unreachable'); addPath(g.from); steps.push({to: t, ...g, det: patch[t].det}); have.add(t); };
  need.forEach(addPath);
  const cellReps = Array.from({length: tiles.length}, (_, c) => { const p = patch[repOf.get(c)]; return {at: repOf.get(c), I: p.I, J: p.J}; });

  const angEq = basisRows(angRows, n + 1), sideEq = basisRows(sideRows, n);
  out[fam.key] = {n, ord, k: fam.k, angles: angEq, sides: sideEq, steps: steps.map(s => [s.to, s.from, s.a, s.a2, s.b, s.b2, s.det]),
    cell: cellReps.map(r => [r.at, r.I, r.J]), lat: [idx(0, 1, 0), idx(0, 0, 1)],
    describe: describe(angEq, sideEq, n)};
  const dofA = n - angEq.length, dofS = n - sideEq.length - 2;
  console.log(fam.key.padEnd(4), `dof ${dofA + dofS - 1}`.padEnd(7), out[fam.key].describe.join('; '));
}

writeFileSync(new URL('./edgesym-data.js', import.meta.url),
  '// Generated by edgesym-derive.mjs: for each type, the linear constraints its tiling\'s edge symmetries impose on\n' +
  '// interior angles (units of pi, vertex i) and side lengths (side i = vertex i -> i+1), and the gluing steps that\n' +
  '// rebuild the cell from the tile: [to, from, a, a2, b, b2, det] = corner b of tile "to" on corner a of "from", ray b->b2 along a->a2.\n' +
  'export const EDGESYM = ' + JSON.stringify(out) + ';\n');

// ---- helpers
function rref(rows, upto = Infinity) {
  const A = rows.map(r => r.slice()); if (!A.length) return A;
  const m = A[0].length; let r = 0;
  for (let c = 0; c < m && r < A.length; c++) {
    let p = r; for (let i = r + 1; i < A.length; i++) if (Math.abs(A[i][c]) > Math.abs(A[p][c])) p = i;
    if (Math.abs(A[p][c]) < 1e-9) continue;
    [A[r], A[p]] = [A[p], A[r]];
    const v = A[r][c]; A[r] = A[r].map(x => x / v);
    for (let i = 0; i < A.length; i++) if (i !== r && Math.abs(A[i][c]) > 0) { const f = A[i][c]; A[i] = A[i].map((x, k) => x - f * A[r][k]); }
    r++;
  }
  return A;
}
// independent rows, reduced, with small-integer coefficients
function basisRows(rows, m) {
  const R = rref(rows).filter(r => r.some(v => Math.abs(v) > 1e-9));
  return R.map(r => {
    for (let d = 1; d <= 12; d++) { const s = r.map(v => v * d); if (s.every(v => Math.abs(v - Math.round(v)) < 1e-7)) return s.map(v => Math.round(v) || 0); }
    throw new Error('non-rational row ' + r);
  });
}
function describe(ang, side, n) {
  const V = 'ABCDEF', S = 'abcdef', term = (c, name) => (c < 0 ? '− ' : '+ ') + (Math.abs(c) === 1 ? '' : Math.abs(c)) + name;
  const lin = (r, names) => r.map((c, i) => c ? term(c, names[i]) : '').filter(Boolean).join(' ').replace(/^\+ /, '').replace(/^− /, '−');
  const A = ang.map(r => lin(r.slice(0, n), [...V].map(v => '∠' + v)) + ' = ' + Math.round(r[n] * 180) + '°');
  const B = side.map(r => {
    const nm = i => S[(i + 1) % n]; // side i runs from vertex i to i+1; named, as in families.js, after its end vertex
    const pos = r.map((c, i) => c > 0 ? (c === 1 ? '' : c) + nm(i) : '').filter(Boolean).join(' + ');
    const neg = r.map((c, i) => c < 0 ? (c === -1 ? '' : -c) + nm(i) : '').filter(Boolean).join(' + ');
    return `${pos} = ${neg || 0}`;
  });
  return [...A, ...B];
}
