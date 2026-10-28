// Periodic tilings by one convex tile: wallpaper group G (p2, pgg, p3, p4, p6) acting on k tile orbits.
// A tiling is a vector x = [lattice params..., (x, y, theta) per orbit]; it is valid when total overlap area is 0
// (the cell area is fixed to k * |G/T| * area(tile), so zero overlap means perfect coverage).

const S3 = Math.sqrt(3);
export const GROUPS = {
  p2:  {order: 2, nlp: 2, name: '2222'},
  pgg: {order: 4, nlp: 1, name: '22×'},
  p3:  {order: 3, nlp: 0, name: '333'},
  p4:  {order: 4, nlp: 0, name: '442'},
  p6:  {order: 6, nlp: 0, name: '632'},
};

export function polyArea(P) {
  let s = 0;
  for (let i = 0; i < P.length; i++) { const a = P[i], b = P[(i + 1) % P.length]; s += a.x * b.y - b.x * a.y; }
  return s / 2;
}
const ccw = P => polyArea(P) < 0 ? P.slice().reverse() : P;

// overlap area of two convex CCW polygons (Sutherland–Hodgman)
export function overlapArea(A, B, dx = 0, dy = 0) {
  let out = A;
  for (let i = 0; i < B.length && out.length; i++) {
    const c1x = B[i].x + dx, c1y = B[i].y + dy, c2x = B[(i + 1) % B.length].x + dx, c2y = B[(i + 1) % B.length].y + dy;
    const ex = c2x - c1x, ey = c2y - c1y;
    const inp = out; out = [];
    for (let j = 0; j < inp.length; j++) {
      const p = inp[j], q = inp[(j + 1) % inp.length];
      const sp = ex * (p.y - c1y) - ey * (p.x - c1x), sq = ex * (q.y - c1y) - ey * (q.x - c1x);
      if (sp >= 0) {
        out.push(p);
        if (sq < 0) { const t = sp / (sp - sq); out.push({x: p.x + t * (q.x - p.x), y: p.y + t * (q.y - p.y)}); }
      } else if (sq >= 0) { const t = sp / (sp - sq); out.push({x: p.x + t * (q.x - p.x), y: p.y + t * (q.y - p.y)}); }
    }
  }
  return out.length < 3 ? 0 : Math.abs(polyArea(out));
}

export function cellArea(cfg, P) { return cfg.k * GROUPS[cfg.group].order * Math.abs(polyArea(P)); }

export function latticeOf(cfg, A, lp) {
  switch (cfg.group) {
    case 'p2': { const p = Math.abs(lp[0]) + 1e-9; return {t1: {x: p, y: 0}, t2: {x: lp[1], y: A / p}}; }
    case 'pgg': { const w = Math.abs(lp[0]) + 1e-9; return {t1: {x: w, y: 0}, t2: {x: 0, y: A / w}}; }
    case 'p4': { const s = Math.sqrt(A); return {t1: {x: s, y: 0}, t2: {x: 0, y: s}}; }
    default: { const s = Math.sqrt(2 * A / S3); return {t1: {x: s, y: 0}, t2: {x: s / 2, y: s * S3 / 2}}; }
  }
}
const rot = a => { const c = Math.cos(a), s = Math.sin(a); return p => ({x: c * p.x - s * p.y, y: s * p.x + c * p.y}); };
function elements(g, lat) {
  switch (g) {
    case 'p2': return [p => p, rot(Math.PI)];
    case 'pgg': { const w = lat.t1.x, h = lat.t2.y;
      return [p => p, p => ({x: -p.x, y: -p.y}), p => ({x: p.x + w / 2, y: -p.y + h / 2}), p => ({x: -p.x + w / 2, y: p.y + h / 2})]; }
    case 'p3': return [0, 1, 2].map(i => rot(i * 2 * Math.PI / 3));
    case 'p4': return [0, 1, 2, 3].map(i => rot(i * Math.PI / 2));
    default: return [0, 1, 2, 3, 4, 5].map(i => rot(i * Math.PI / 3));
  }
}

// cell contents: lattice plus every tile image (orbit index, CCW polygon)
export function buildCell(cfg, P, x) {
  const A = cellArea(cfg, P), nl = GROUPS[cfg.group].nlp;
  const lat = latticeOf(cfg, A, x.slice(0, nl));
  const els = elements(cfg.group, lat), tiles = [];
  for (let j = 0; j < cfg.k; j++) {
    let px = x[nl + 3 * j], py = x[nl + 3 * j + 1];
    const th = x[nl + 3 * j + 2];
    { // reduce the base position modulo the lattice (images shift by lattice vectors too)
      const det = lat.t1.x * lat.t2.y - lat.t1.y * lat.t2.x;
      const fa = Math.round((px * lat.t2.y - py * lat.t2.x) / det - 0.137), fb = Math.round((lat.t1.x * py - lat.t1.y * px) / det - 0.211); // offsets keep exact half-lattice positions off the rounding boundary
      px -= fa * lat.t1.x + fb * lat.t2.x; py -= fa * lat.t1.y + fb * lat.t2.y;
    }
    const R = rot(th), f = cfg.flips[j] ? -1 : 1;
    const base = P.map(v => { const r = R({x: v.x, y: f * v.y}); return {x: r.x + px, y: r.y + py}; });
    for (const g of els) tiles.push({orbit: j, poly: ccw(base.map(g))});
  }
  return {lat, tiles};
}

function reduceBasis(a, b) { // Gauss reduction: short, nearly orthogonal basis of the same lattice
  for (let it = 0; it < 50; it++) {
    if (a.x * a.x + a.y * a.y > b.x * b.x + b.y * b.y) [a, b] = [b, a];
    const m = Math.round((a.x * b.x + a.y * b.y) / (a.x * a.x + a.y * a.y));
    if (!m) break;
    b = {x: b.x - m * a.x, y: b.y - m * a.y};
  }
  return {a, b};
}
// ear-clipping triangulation of a simple CCW polygon (for non-convex tiles)
export function triangulate(poly) {
  const idx = poly.map((_, i) => i), tris = [];
  const cr = (a, b, c) => (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  while (idx.length > 3) {
    let clipped = false;
    for (let i = 0; i < idx.length && !clipped; i++) {
      const ia = idx[(i + idx.length - 1) % idx.length], ib = idx[i], ic = idx[(i + 1) % idx.length];
      const a = poly[ia], b = poly[ib], c = poly[ic];
      if (cr(a, b, c) <= 1e-12) continue;
      if (idx.some(j => j !== ia && j !== ib && j !== ic && cr(a, b, poly[j]) >= -1e-12 && cr(b, c, poly[j]) >= -1e-12 && cr(c, a, poly[j]) >= -1e-12)) continue;
      tris.push([a, b, c]); idx.splice(i, 1); clipped = true;
    }
    if (!clipped) break;
  }
  if (idx.length === 3) tris.push(idx.map(i => poly[i]));
  return tris;
}

export function energy(cfg, P, x, nonconvex = false) {
  const {lat, tiles} = buildCell(cfg, P, x);
  const n = tiles.length, info = tiles.map(t => {
    let cx = 0, cy = 0; for (const v of t.poly) { cx += v.x; cy += v.y; } cx /= t.poly.length; cy /= t.poly.length;
    let r = 0; for (const v of t.poly) r = Math.max(r, Math.hypot(v.x - cx, v.y - cy));
    return {cx, cy, r};
  });
  let e = 0;
  const tris = [];
  const R = 3, {a: r1, b: r2} = reduceBasis(lat.t1, lat.t2);
  for (let i = -R; i <= R; i++) for (let j = -R; j <= R; j++) {
    if (i < 0 || (i === 0 && j < 0)) continue;
    const sx = i * r1.x + j * r2.x, sy = i * r1.y + j * r2.y;
    for (let u = 0; u < n; u++) for (let v = 0; v < n; v++) {
      if (i === 0 && j === 0 && v <= u) continue;
      if (Math.hypot(info[v].cx + sx - info[u].cx, info[v].cy + sy - info[u].cy) >= info[u].r + info[v].r) continue;
      if (!nonconvex) { e += overlapArea(tiles[u].poly, tiles[v].poly, sx, sy); continue; }
      tris[u] = tris[u] || triangulate(tiles[u].poly); tris[v] = tris[v] || triangulate(tiles[v].poly);
      for (const ta of tris[u]) for (const tb of tris[v]) e += overlapArea(ta, tb, sx, sy);
    }
  }
  return e / Math.abs(polyArea(P));
}

export function nelderMead(f, x0, step, maxEval, tol = 1e-12) {
  const n = x0.length;
  let pts = [x0.slice()];
  for (let i = 0; i < n; i++) { const p = x0.slice(); p[i] += step[i]; pts.push(p); }
  let vals = pts.map(f), evals = n + 1;
  while (evals < maxEval) {
    const order = vals.map((v, i) => i).sort((a, b) => vals[a] - vals[b]);
    pts = order.map(i => pts[i]); vals = order.map(i => vals[i]);
    if (vals[0] < tol) break;
    const c = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) c[k] += pts[i][k] / n;
    const at = (t) => c.map((v, k) => v + t * (pts[n][k] - v));
    const xr = at(-1), fr = f(xr); evals++;
    if (fr < vals[0]) {
      const xe = at(-2), fe = f(xe); evals++;
      if (fe < fr) { pts[n] = xe; vals[n] = fe; } else { pts[n] = xr; vals[n] = fr; }
    } else if (fr < vals[n - 1]) { pts[n] = xr; vals[n] = fr; }
    else {
      const xc = fr < vals[n] ? at(-0.5) : at(0.5), fc = f(xc); evals++;
      if (fc < Math.min(fr, vals[n])) { pts[n] = xc; vals[n] = fc; }
      else for (let i = 1; i <= n; i++) { pts[i] = pts[i].map((v, k) => pts[0][k] + 0.5 * (v - pts[0][k])); vals[i] = f(pts[i]); evals += 1; }
    }
  }
  let b = 0; for (let i = 1; i <= n; i++) if (vals[i] < vals[b]) b = i;
  return {x: pts[b], e: vals[b]};
}

// refine x toward zero overlap, starting from a nearby solution
export function refine(cfg, P, x, scale, maxEval = 3000) {
  const f = z => energy(cfg, P, z);
  let best = {x, e: f(x)};
  for (let r = 0; r < 4 && best.e > 1e-11; r++) {
    const s = best.x.map(() => scale * 0.05 / (1 + 3 * r));
    const res = nelderMead(f, best.x, s, maxEval, 1e-12);
    if (res.e < best.e) best = res;
  }
  return best;
}

export function randomStart(cfg, P, rnd) {
  const A = cellArea(cfg, P), s = Math.sqrt(A), x = [];
  if (cfg.group === 'p2') { const p = s * Math.exp((rnd() - .5) * 1.0); x.push(p, (rnd() - .5) * p); }
  if (cfg.group === 'pgg') x.push(s * Math.exp((rnd() - .5) * 1.2));
  for (let j = 0; j < cfg.k; j++) x.push((rnd() - .5) * s, (rnd() - .5) * s, rnd() * 2 * Math.PI);
  return x;
}

export function solveSym(M, b) { // Gaussian elimination, M is n x n
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

// ---- contacts: vertex-on-edge incidences of a (near-)perfect tiling, kept as equations for continuation
function distPtSeg(p, a, b) {
  const ex = b.x - a.x, ey = b.y - a.y, t = Math.max(0, Math.min(1, ((p.x - a.x) * ex + (p.y - a.y) * ey) / (ex * ex + ey * ey)));
  return Math.hypot(p.x - a.x - t * ex, p.y - a.y - t * ey);
}
export function findContacts(cfg, P, x, tol) {
  const {lat, tiles} = buildCell(cfg, P, x), out = [];
  for (let i = -3; i <= 3; i++) for (let j = -3; j <= 3; j++) {
    const sx = i * lat.t1.x + j * lat.t2.x, sy = i * lat.t1.y + j * lat.t2.y;
    tiles.forEach((tu, u) => tiles.forEach((tv, v) => {
      if (i === 0 && j === 0 && u === v) return;
      const m = tv.poly.length;
      tu.poly.forEach((p, a) => {
        for (let b = 0; b < m; b++) {
          const c1 = {x: tv.poly[b].x + sx, y: tv.poly[b].y + sy}, c2 = {x: tv.poly[(b + 1) % m].x + sx, y: tv.poly[(b + 1) % m].y + sy};
          if (distPtSeg(p, c1, c2) < tol) out.push([u, a, v, b, i, j]);
        }
      });
    }));
  }
  return out;
}
export function contactResiduals(cfg, P, x, contacts) {
  const {lat, tiles} = buildCell(cfg, P, x);
  return contacts.map(([u, a, v, b, i, j]) => {
    const sx = i * lat.t1.x + j * lat.t2.x, sy = i * lat.t1.y + j * lat.t2.y, q = tiles[v].poly, m = q.length;
    const c1 = q[b], c2 = q[(b + 1) % m], p = tiles[u].poly[a];
    const ex = c2.x - c1.x, ey = c2.y - c1.y;
    return (ex * (p.y - c1.y - sy) - ey * (p.x - c1.x - sx)) / Math.hypot(ex, ey);
  });
}
// Levenberg–Marquardt on the contact equations; returns {x, r} (r = max residual)
export function solveContacts(cfg, P, x, contacts, iters = 25) {
  x = x.slice(); const n = x.length, h = 1e-7;
  const maxAbs = r => r.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
  let r = contactResiduals(cfg, P, x, contacts), lam = 1e-9;
  for (let it = 0; it < iters && maxAbs(r) > 1e-12 && lam < 1e6; it++) {
    const J = Array.from({length: n}, (_, k) => { const y = x.slice(); y[k] += h; return contactResiduals(cfg, P, y, contacts).map((v, q) => (v - r[q]) / h); });
    const M = J.map((jr, a) => J.map((jc, b) => jr.reduce((s, v, q) => s + v * jc[q], 0) + (a === b ? lam : 0)));
    const g = J.map(jr => -jr.reduce((s, v, q) => s + v * r[q], 0));
    const dx = solveSym(M, g), y = x.map((v, k) => v + dx[k]), r2 = contactResiduals(cfg, P, y, contacts);
    if (r2.reduce((a, v) => a + v * v, 0) < r.reduce((a, v) => a + v * v, 0)) { x = y; r = r2; lam = Math.max(lam / 3, 1e-14); } else lam *= 10;
  }
  return {x, r: maxAbs(r)};
}
// try to turn a near-miss into an exact tiling
export function polish(cfg, P, x) {
  for (const tol of [1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3]) {
    const c = findContacts(cfg, P, x, tol);
    if (!c.length) continue;
    const s = solveContacts(cfg, P, x, c);
    const e = energy(cfg, P, s.x);
    if (e < 1e-9) return {x: s.x, e};
  }
  return null;
}

// start with each tile touching an earlier one (vertex on vertex, edges aligned)
export function structuredStart(cfg, P, rnd) {
  const x = randomStart(cfg, P, rnd), nl = GROUPS[cfg.group].nlp, n = P.length;
  for (let j = 0; j < cfg.k; j++) {
    const f = cfg.flips[j] ? -1 : 1, own = P.map(v => ({x: v.x, y: f * v.y}));
    let th, px, py;
    if (j === 0) {
      th = rnd() * 2 * Math.PI;
      const a = Math.floor(rnd() * 2 * n), b = a < n ? own[a] : {x: (own[a - n].x + own[(a - n + 1) % n].x) / 2, y: (own[a - n].y + own[(a - n + 1) % n].y) / 2};
      const r = rot(th)(b); px = -r.x; py = -r.y;
    } else {
      const part = buildCell({...cfg, k: j}, P, x.slice(0, nl + 3 * j)).tiles;
      const Q = part[Math.floor(rnd() * part.length)].poly, qa = Math.floor(rnd() * n), a = Math.floor(rnd() * n);
      const qd = rnd() < .5 ? 1 : -1, od = rnd() < .5 ? 1 : -1;
      const d = {x: Q[(qa + qd + n) % n].x - Q[qa].x, y: Q[(qa + qd + n) % n].y - Q[qa].y};
      const u = {x: own[(a + od + n) % n].x - own[a].x, y: own[(a + od + n) % n].y - own[a].y};
      th = Math.atan2(d.y, d.x) + Math.PI - Math.atan2(u.y, u.x);
      const r = rot(th)(own[a]); px = Q[qa].x - r.x; py = Q[qa].y - r.y;
    }
    x[nl + 3 * j] = px; x[nl + 3 * j + 1] = py; x[nl + 3 * j + 2] = th;
  }
  return x;
}

// sample many random configurations, refine the most promising; returns {x, e} (e < 1e-9 means a perfect tiling)
export function findTiling(cfg, P, rnd, samples, keep = 30) {
  const A = cellArea(cfg, P), scale = Math.sqrt(A);
  const f = z => energy(cfg, P, z);
  const pool = [];
  for (let t = 0; t < samples; t++) { const x = (t & 1) ? structuredStart(cfg, P, rnd) : randomStart(cfg, P, rnd); pool.push({x, e: f(x)}); }
  pool.sort((a, b) => a.e - b.e);
  let bestAll = null;
  for (const c of pool.slice(0, keep)) {
    let res = c;
    for (let r = 0; r < 8 && res.e > 1e-10; r++) {
      const r2 = nelderMead(f, res.x, res.x.map(() => scale * 0.12 / (1 + r)), 1200, 1e-12);
      if (r2.e < res.e) res = r2;
    }
    if (!bestAll || res.e < bestAll.e) bestAll = res;
    if (res.e < 1e-9) return res;
    if (res.e < 1e-2) { const p = polish(cfg, P, res.x); if (p) return p; }
  }
  return bestAll;
}

export function mulberry(seed) {
  return () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
}
