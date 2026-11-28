// Edge-symmetry parametrization (experiment). Each type's tiling pairs edges and corners by its symmetries; that
// gives linear constraints on interior angles and on side lengths (edgesym-data.js). We parametrize that space
// directly: angles = g0 + Na·u (affine), sides = Bs·w made to close up by a linear projection (closure is linear in the
// sides once the angles are fixed). The tile is then laid down by turtle, and the cell is built by gluing tiles corner
// to corner along shared edge rays — no projection onto a constraint manifold, no tiling solver.
import {EDGESYM} from './edgesym-data.js';
import {overlapArea, triangulate, reduceBasis} from './tiler.js';
export {EDGESYM};

const PI = Math.PI;
const dot = (a, b) => a.reduce((s, v, i) => s + v * b[i], 0);
// orthonormal basis of the row space of A (rows), via Gram–Schmidt
function rowBasis(A) {
  const B = [];
  for (const r of A) {
    let v = r.slice(); for (const b of B) { const c = dot(v, b); v = v.map((x, i) => x - c * b[i]); }
    const l = Math.hypot(...v); if (l > 1e-9) B.push(v.map(x => x / l));
  }
  return B;
}
function nullBasis(A, m) {
  const Rb = rowBasis(A), out = [];
  for (let k = 0; k < m; k++) {
    let v = Array.from({length: m}, (_, i) => +(i === k));
    for (const b of [...Rb, ...out]) { const c = dot(v, b); v = v.map((x, i) => x - c * b[i]); }
    const l = Math.hypot(...v); if (l > 1e-9) out.push(v.map(x => x / l));
  }
  return out;
}
const anglesOf = P => { const n = P.length / 2; return Array.from({length: n}, (_, i) => {
  const p = (i + n - 1) % n, q = (i + 1) % n;
  const ax = P[2*i] - P[2*p], ay = P[2*i+1] - P[2*p+1], bx = P[2*q] - P[2*i], by = P[2*q+1] - P[2*i+1];
  return PI - Math.atan2(ax * by - ay * bx, ax * bx + ay * by);
}); };
const sidesOf = P => { const n = P.length / 2; return Array.from({length: n}, (_, i) => { const j = (i + 1) % n; return Math.hypot(P[2*j] - P[2*i], P[2*j+1] - P[2*i+1]); }); };

export function makeES(key, P0) {
  const D = EDGESYM[key], n = D.n;
  const Aa = D.angles.map(r => r.slice(0, n)), ba = D.angles.map(r => r[n] * PI);
  // particular angle solution nearest the stored tile's angles
  let g0 = anglesOf(P0);
  { // min-norm correction onto the constraints, via the normal equations
    const res = Aa.map((r, i) => ba[i] - dot(r, g0)), G = Aa.map(r => Aa.map(q => dot(r, q)));
    const y = gauss(G, res); g0 = g0.map((v, k) => v + Aa.reduce((s, r, i) => s + r[k] * y[i], 0)); }
  const Na = nullBasis(Aa, n), Bs = nullBasis(D.sides, n), da = Na.length, ds = Bs.length;
  const dof = da + ds - 2 - 1;

  const angles = u => g0.map((v, i) => v + Na.reduce((s, b, k) => s + b[i] * u[k], 0));
  // theta = [u (da), w (ds), x0, y0, heading of side 0]
  const heads = (u, rot) => { const g = angles(u), h = [rot]; for (let i = 1; i < n; i++) h.push(h[i - 1] + PI - g[i]); return {g, h}; };
  // closure vectors of the side basis: z_k = sum_i Bs[k][i] e^{i h_i}
  const zs = h => Bs.map(b => [dot(b, h.map(Math.cos)), dot(b, h.map(Math.sin))]);
  function shape(th) {
    let u = th.slice(0, da), w = th.slice(da, da + ds); const [x0, y0, rot] = th.slice(da + ds);
    if (ds === 2) {
      // only two side directions: closure w1 z1 + w2 z2 = 0 needs z1 parallel to z2 — one real equation on the angles
      const f = v => { const [a, b] = zs(heads(v, 0).h); return (a[0] * b[1] - a[1] * b[0]) / (Math.hypot(...a) * Math.hypot(...b)); };
      for (let it = 0; it < 30; it++) {
        const f0 = f(u); if (Math.abs(f0) < 1e-15) break;
        const gr = u.map((_, k) => { const v = u.slice(); v[k] += 1e-7; return (f(v) - f0) / 1e-7; }), gg = dot(gr, gr);
        if (!gg) break; u = u.map((v, k) => v - f0 * gr[k] / gg);
      }
      const [a, b] = zs(heads(u, 0).h), lam = dot(a, b) / dot(b, b), r = Math.hypot(...w) || 1;
      const w2 = [1, -lam];
      const sum = Array.from({length: n}, (_, i) => Bs[0][i] * w2[0] + Bs[1][i] * w2[1]).reduce((t, v) => t + v, 0);
      const k = r / Math.hypot(...w2) * Math.sign(sum || 1); w = w2.map(v => v * k);
    }
    const {g, h} = heads(u, rot);
    const K = ds === 2 ? [] : rowBasis([h.map(Math.cos), h.map(Math.sin)].map(c => Bs.map(b => dot(c, b))));
    let wc = w.slice(); for (const k of K) { const c = dot(wc, k); wc = wc.map((x, i) => x - c * k[i]); }
    const s = Array.from({length: n}, (_, i) => Bs.reduce((t, b, k) => t + b[i] * wc[k], 0));
    const P = [x0, y0];
    for (let i = 0; i < n - 1; i++) P.push(P[2*i] + s[i] * Math.cos(h[i]), P[2*i+1] + s[i] * Math.sin(h[i]));
    return {P, g, s};
  }
  function thetaOf(P) {
    const g = anglesOf(P), s = sidesOf(P);
    return [...Na.map(b => dot(b, g.map((v, i) => v - g0[i]))), ...Bs.map(b => dot(b, s)), P[0], P[1], Math.atan2(P[3] - P[1], P[2] - P[0])];
  }
  // exact member of the space nearest P
  const snap = P => shape(thetaOf(P)).P;

  // move vertex i toward q (tile frame), disturbing the other vertices as little as possible
  function drag(P, i, q, mu = 0.05) {
    let th = thetaOf(P); const m = th.length;
    const res = t => { const Q = shape(t).P, r = [Q[2*i] - q[0], Q[2*i+1] - q[1]];
      for (let j = 0; j < n; j++) if (j !== i) r.push(mu * (Q[2*j] - P[2*j]), mu * (Q[2*j+1] - P[2*j+1]));
      return r; };
    let r = res(th), lam = 1e-6;
    for (let it = 0; it < 12; it++) {
      const J = Array.from({length: m}, (_, k) => { const t = th.slice(); t[k] += 1e-7; return res(t).map((v, j) => (v - r[j]) / 1e-7); });
      const M = J.map((a, x) => J.map((b, y) => dot(a, b) + (x === y ? lam : 0))), gr = J.map(a => -dot(a, r));
      const t2 = th.map((v, k) => v + gauss(M, gr)[k]), r2 = res(t2);
      if (dot(r2, r2) < dot(r, r)) { th = t2; r = r2; lam = Math.max(lam / 3, 1e-12); } else lam *= 10;
      if (dot(r, r) < 1e-24) break;
    }
    return shape(th);
  }

  // glue the cell together from the tile
  function build(P) {
    const V = Array.from({length: n}, (_, i) => ({x: P[2*i], y: P[2*i+1]}));
    const placed = new Map();
    const id = {L: [1, 0, 0, 1], t: {x: 0, y: 0}}, ap = (M, v) => ({x: M.L[0] * v.x + M.L[1] * v.y + M.t.x, y: M.L[2] * v.x + M.L[3] * v.y + M.t.y});
    const start = D.steps.length ? D.steps[0][1] : D.cell[0][0];
    placed.set(start, {M: id, W: V, det: 1});
    for (const [to, from, a, a2, b, b2, det] of D.steps) {
      const A = placed.get(from), pa = A.W[a], va = Math.atan2(A.W[a2].y - pa.y, A.W[a2].x - pa.x);
      const vb = Math.atan2(V[b2].y - V[b].y, V[b2].x - V[b].x);
      const al = det > 0 ? va - vb : va + vb, c = Math.cos(al), s = Math.sin(al);
      const L = det > 0 ? [c, -s, s, c] : [c, s, s, -c];
      const M = {L, t: {x: 0, y: 0}}, pb = ap(M, V[b]); M.t = {x: pa.x - pb.x, y: pa.y - pb.y};
      placed.set(to, {M, W: V.map(v => ap(M, v)), det});
    }
    const W0 = placed.get(start).W[0];
    const t1 = {x: placed.get(D.lat[0]).W[0].x - W0.x, y: placed.get(D.lat[0]).W[0].y - W0.y};
    const t2 = {x: placed.get(D.lat[1]).W[0].x - W0.x, y: placed.get(D.lat[1]).W[0].y - W0.y};
    const tiles = [], els = [];
    D.cell.forEach(([at, I, J], c) => {
      const T = placed.get(at), sx = -(I * t1.x + J * t2.x), sy = -(I * t1.y + J * t2.y);
      const W = T.W.map(v => ({x: v.x + sx, y: v.y + sy}));
      tiles.push({orbit: Math.floor(c / D.ord), poly: T.det > 0 ? W : W.slice().reverse()});
      if (c < D.ord) els.push({L: T.M.L, t: {x: T.M.t.x + sx, y: T.M.t.y + sy}});
    });
    return {lat: {t1, t2}, tiles, els};
  }
  return {key, dof, shape, thetaOf, snap, drag, build, describe: D.describe};
}

// total overlap area of a cell's tiles with each other and their lattice translates, relative to one tile's area
export function cellOverlap(cell) {
  const {tiles, lat} = cell, info = tiles.map(t => {
    let cx = 0, cy = 0; for (const v of t.poly) { cx += v.x; cy += v.y; } cx /= t.poly.length; cy /= t.poly.length;
    let r = 0; for (const v of t.poly) r = Math.max(r, Math.hypot(v.x - cx, v.y - cy));
    return {cx, cy, r, tris: triangulate(t.poly)};
  });
  let area = 0; for (const tr of info[0].tris) area += Math.abs((tr[1].x - tr[0].x) * (tr[2].y - tr[0].y) - (tr[2].x - tr[0].x) * (tr[1].y - tr[0].y)) / 2;
  const {a: r1, b: r2} = reduceBasis(lat.t1, lat.t2), N = tiles.length;
  let e = 0;
  for (let i = -3; i <= 3; i++) for (let j = -3; j <= 3; j++) {
    if (i < 0 || (i === 0 && j < 0)) continue;
    const sx = i * r1.x + j * r2.x, sy = i * r1.y + j * r2.y;
    for (let u = 0; u < N; u++) for (let v = 0; v < N; v++) {
      if (i === 0 && j === 0 && v <= u) continue;
      if (Math.hypot(info[v].cx + sx - info[u].cx, info[v].cy + sy - info[u].cy) >= info[u].r + info[v].r) continue;
      for (const ta of info[u].tris) for (const tb of info[v].tris) e += overlapArea(ta, tb, sx, sy);
    }
  }
  // a gap shows up as the cell being too large: compare cell area with the tiles' total area
  const det = Math.abs(lat.t1.x * lat.t2.y - lat.t1.y * lat.t2.x);
  return (e + Math.abs(det - N * area)) / area;
}

function gauss(M, b) {
  const n = b.length, A = M.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
    [A[c], A[p]] = [A[p], A[c]];
    const v = A[c][c] || 1e-300;
    for (let r = c + 1; r < n; r++) { const f = A[r][c] / v; for (let k = c; k <= n; k++) A[r][k] -= f * A[c][k]; }
  }
  const x = new Array(n).fill(0);
  for (let r = n - 1; r >= 0; r--) { let s = A[r][n]; for (let k = r + 1; k < n; k++) s -= A[r][k] * x[k]; x[r] = s / (A[r][r] || 1e-300); }
  return x;
}
