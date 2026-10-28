// Reinhardt-style families of convex monotiles. Vertices A,B,C,...
// Vertices run counterclockwise; side a is E->A, b: A->B, c: B->C, ... (so A is opposite d, B opposite e, ...).
// A family is a set of equations on the vertex coordinates (flat array [x0,y0,x1,y1,...]).
const d = Math.PI / 180, PI = Math.PI;

const sideLen = (P, i) => { const n = P.length / 2, j = (i + 1) % n; return Math.hypot(P[2*j] - P[2*i], P[2*j+1] - P[2*i+1]); };
// interior angle at vertex i of a counterclockwise polygon
const angleAt = (P, i) => {
  const n = P.length / 2, p = (i + n - 1) % n, q = (i + 1) % n;
  const ax = P[2*i] - P[2*p], ay = P[2*i+1] - P[2*p+1], bx = P[2*q] - P[2*i], by = P[2*q+1] - P[2*i+1];
  return PI - Math.atan2(ax * by - ay * bx, ax * bx + ay * by);
};
const SH = {v: -1};
const mk = f => P => { const n = P.length / 2, s = i => sideLen(P, ((i + SH.v) % n + n) % n), g = i => angleAt(P, i); return f(s, g); };

export const FAMILIES = [
  {key: 'tri', n: 3, title: 'Triangle', label: 'any triangle', group: 'p2', k: 1, flips: [0], cons: mk(() => [])},
  {key: 'quad', n: 4, title: 'Quadrilateral', label: 'any quadrilateral', group: 'p2', k: 1, flips: [0], cons: mk(() => [])},

  {key: 'p1', n: 5, label: 'Type 1', group: 'p2', k: 1, flips: [0], cons: mk((s, g) => [g(1) + g(2) - PI])},
  {key: 'p2', n: 5, label: 'Type 2', group: 'pgg', k: 1, flips: [0], cons: mk((s, g) => [s(2) - s(4), g(1) + g(3) - PI])},
  {key: 'p3', n: 5, label: 'Type 3', group: 'p3', k: 1, flips: [0], cons: mk((s, g) => [s(0) - s(1), s(3) - s(2) - s(4), g(0) - 120*d, g(2) - 120*d, g(3) - 120*d])},
  {key: 'p4', n: 5, label: 'Type 4', group: 'p4', k: 1, flips: [0], cons: mk((s, g) => [s(1) - s(2), s(3) - s(4), g(1) - 90*d, g(3) - 90*d])},
  {key: 'p5', n: 5, label: 'Type 5', group: 'p6', k: 1, flips: [0], cons: mk((s, g) => [s(0) - s(1), s(3) - s(4), g(0) - 60*d, g(3) - 120*d])},
  {key: 'p6', n: 5, label: 'Type 6', group: 'p2', k: 2, flips: [0, 0], cons: mk((s, g) => [s(0) - s(3), s(0) - s(4), s(1) - s(2), g(1) + g(3) - PI, 2*g(1) - g(4)])},
  {key: 'p7', n: 5, label: 'Type 7', group: 'pgg', k: 2, flips: [0, 0], cons: mk((s, g) => [s(1) - s(2), s(1) - s(3), s(1) - s(4), g(1) + 2*g(4) - 2*PI, 2*g(2) + g(3) - 2*PI])},
  {key: 'p8', n: 5, label: 'Type 8', group: 'pgg', k: 2, flips: [0, 0], cons: mk((s, g) => [s(1) - s(2), s(1) - s(3), s(1) - s(4), 2*g(1) + g(2) - 2*PI, g(3) + 2*g(4) - 2*PI])},
  {key: 'p9', n: 5, label: 'Type 9', group: 'pgg', k: 2, flips: [0, 0], cons: mk((s, g) => [s(1) - s(2), s(1) - s(3), s(1) - s(4), 2*g(0) + g(2) - 2*PI, g(3) + 2*g(4) - 2*PI])},
  {key: 'p10', n: 5, label: 'Type 10', group: 'p2', k: 3, flips: [0, 0, 0], cons: mk((s, g) => [s(0) - s(1), s(0) - s(2) - s(4), g(0) - 90*d, g(1) + g(4) - PI, g(1) + 2*g(2) - 2*PI])},
  {key: 'p11', n: 5, label: 'Type 11', group: 'pgg', k: 2, flips: [0, 0], cons: mk((s, g) => [2*s(0) + s(2) - s(3), s(3) - s(4), g(0) - 90*d, g(2) + g(4) - PI, 2*g(1) + g(2) - 2*PI])},
  {key: 'p12', n: 5, label: 'Type 12', group: 'pgg', k: 2, flips: [0, 0], cons: mk((s, g) => [2*s(0) - s(3), s(3) - s(2) - s(4), g(0) - 90*d, g(2) + g(4) - PI, 2*g(1) + g(2) - 2*PI])},
  {key: 'p13', n: 5, label: 'Type 13', group: 'pgg', k: 2, flips: [0, 0], cons: mk((s, g) => [s(3) - 2*s(0), s(3) - 2*s(4), g(1) - 90*d, g(4) - 90*d, 2*g(0) + g(3) - 2*PI])},
  {key: 'p14', n: 5, label: 'Type 14', group: 'p2', k: 3, flips: [0, 0, 0], cons: mk((s, g) => [2*s(0) - 2*s(2), s(3) - 2*s(0), s(4) - 2*s(0), g(0) - 90*d, 2*g(1) + g(2) - 2*PI, g(2) + g(4) - PI])},
  {key: 'p15', n: 5, label: 'Type 15', group: 'pgg', k: 3, flips: [0, 0, 0], cons: mk((s, g) => [s(0) - s(2), s(0) - s(4), s(1) - 2*s(0), g(0) - 150*d, g(1) - 60*d, g(2) - 135*d, g(3) - 105*d])},

  {key: 'h1', n: 6, label: 'Type 1', group: 'p2', k: 1, flips: [0], cons: mk((s, g) => [g(0) + g(1) + g(2) - 2*PI, s(0) - s(3)])},
  {key: 'h2', n: 6, label: 'Type 2', group: 'pgg', k: 1, flips: [0], cons: mk((s, g) => [g(0) + g(1) + g(3) - 2*PI, s(0) - s(3), s(2) - s(4)])},
  {key: 'h3', n: 6, label: 'Type 3', group: 'p3', k: 1, flips: [0], cons: mk((s, g) => [g(0) - 120*d, g(2) - 120*d, g(4) - 120*d, s(0) - s(1), s(2) - s(3), s(4) - s(5)])},
];
FAMILIES.forEach(f => { f.title = f.title || {5: 'Pentagon', 6: 'Hexagon'}[f.n]; });

export const setSideShift = v => { SH.v = v; };
export const orbifold = {p2: '2222', pgg: '22×', p3: '333', p4: '442', p6: '632'};

import {solveSym} from './tiler.js';
// ---- small dense linear algebra
export function jacobian(fn, P) {
  const c0 = fn(P), J = c0.map(() => new Array(P.length));
  for (let k = 0; k < P.length; k++) {
    const Q = P.slice(); Q[k] += 1e-6; const c1 = fn(Q);
    for (let r = 0; r < c0.length; r++) J[r][k] = (c1[r] - c0[r]) / 1e-6;
  }
  return {c: c0, J};
}
// Gauss-Newton projection onto the constraint set (min-norm steps)
export function project(fam, P, iters = 30) {
  P = P.slice();
  for (let it = 0; it < iters; it++) {
    const {c, J} = jacobian(fam.cons, P);
    if (!c.length || Math.max(...c.map(Math.abs)) < 1e-12) break;
    const m = c.length, G = J.map((r, i) => J.map((q, j) => r.reduce((s, v, k) => s + v * q[k], 0) + (i === j ? 1e-12 : 0)));
    const y = solveSym(G, c.map(v => -v));
    for (let k = 0; k < P.length; k++) for (let r = 0; r < m; r++) P[k] += J[r][k] * y[r];
  }
  return P;
}
export const residual = (fam, P) => Math.max(0, ...fam.cons(P).map(Math.abs));

// pull vertex i toward target q while staying on the family
export function dragShape(fam, P, i, q) {
  let Q = P.slice();
  for (let it = 0; it < 8; it++) {
    const {c, J} = jacobian(fam.cons, Q), n = Q.length, W = 100;
    const rows = [], rhs = [];
    c.forEach((cv, r) => { rows.push(J[r].map(v => v * W)); rhs.push(-cv * W); });
    for (let a = 0; a < 2; a++) { const r = new Array(n).fill(0); r[2*i+a] = 1; rows.push(r); rhs.push(q[a] - Q[2*i+a]); }
    const M = Array.from({length: n}, (_, a) => Array.from({length: n}, (_, b) => rows.reduce((s, r) => s + r[a] * r[b], 0) + (a === b ? 1e-6 : 0)));
    const v = Array.from({length: n}, (_, a) => rows.reduce((s, r, ri) => s + r[a] * rhs[ri], 0));
    const dx = solveSym(M, v);
    Q = Q.map((x, k) => x + dx[k]);
  }
  return project(fam, Q);
}

export function convexMargin(P) { // min sin(turn)-like margin, <=0 if not strictly convex clockwise; also min edge ratio
  const n = P.length / 2; let m = 1e9, per = 0, tot = 0;
  for (let i = 0; i < n; i++) per += sideLen(P, i);
  for (let i = 0; i < n; i++) {
    const p = (i + n - 1) % n, q = (i + 1) % n;
    const ax = P[2*i] - P[2*p], ay = P[2*i+1] - P[2*p+1], bx = P[2*q] - P[2*i], by = P[2*q+1] - P[2*i+1];
    const la = Math.hypot(ax, ay), lb = Math.hypot(bx, by);
    if (la < 0.06 * per / n || lb < 0.06 * per / n) return -1;
    const cr = ax * by - ay * bx; // counterclockwise => positive
    m = Math.min(m, cr / (la * lb));
    tot += Math.atan2(cr, ax * bx + ay * by);
  }
  if (Math.abs(tot - 2 * PI) > 1e-6) return -1;
  return m;
}
export const interiorAngles = P => Array.from({length: P.length / 2}, (_, i) => angleAt(P, i));

function regular(n, rnd, jitter) {
  const P = [];
  for (let i = 0; i < n; i++) { const a = 2 * PI * i / n, r = 1 + (rnd() - .5) * jitter; P.push(r * Math.cos(a) + (rnd() - .5) * jitter * .6, r * Math.sin(a) + (rnd() - .5) * jitter * .6); }
  return P;
}
// a nicely proportioned member of the family
export function startShape(fam, rnd) {
  if (fam.n === 3) return [1.35, 0.35, 0.45, 1.0, 0, 0];
  if (fam.n === 4) return [1.15, 0.1, 1.0, 0.95, 0.1, 0.8, 0, 0];
  let best = null, bs = -1;
  for (let t = 0; t < 4000; t++) {
    let P = project(fam, regular(fam.n, rnd, 0.9 + rnd() * 0.8), 60);
    if (residual(fam, P) > 1e-9) continue;
    const m = convexMargin(P); if (m <= 0.03) continue;
    const per = Array.from({length: fam.n}, (_, i) => sideLen(P, i)); const ratio = Math.min(...per) / Math.max(...per);
    const score = m * ratio;
    if (score > bs) { bs = score; best = P; }
  }
  return best;
}
export function polyOf(P) { const o = []; for (let i = 0; i < P.length; i += 2) o.push({x: P[i], y: P[i+1]}); return o; }

// number of genuine shape degrees of freedom (beyond translation, rotation, scale)
export function shapeDof(fam, P) {
  const {J} = jacobian(fam.cons, P), n = P.length;
  const A = J.map(r => r.slice()); let rank = 0, big = 0;
  A.forEach(r => r.forEach(v => { big = Math.max(big, Math.abs(v)); }));
  for (let c = 0; c < n && rank < A.length; c++) {
    let p = rank; for (let r = rank + 1; r < A.length; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
    if (!A.length || Math.abs(A[p][c]) < 1e-5 * (big || 1)) continue;
    [A[rank], A[p]] = [A[p], A[rank]];
    for (let r = rank + 1; r < A.length; r++) { const f = A[r][c] / A[rank][c]; for (let k = c; k < n; k++) A[r][k] -= f * A[rank][k]; }
    rank++;
  }
  return n - rank - 4;
}

// a usable tile: simple polygon (convex or not), no tiny edges, no needle-sharp corners
export function tileOK(P) {
  const n = P.length / 2, V = i => ({x: P[2 * ((i + n) % n)], y: P[2 * ((i + n) % n) + 1]});
  let area = 0, per = 0;
  for (let i = 0; i < n; i++) { const a = V(i), b = V(i + 1); area += a.x * b.y - b.x * a.y; per += Math.hypot(b.x - a.x, b.y - a.y); }
  if (area <= 0) return false;
  for (let i = 0; i < n; i++) {
    const a = V(i), b = V(i + 1);
    if (Math.hypot(b.x - a.x, b.y - a.y) < 0.06 * per / n) return false;
    const g = angleAt(P, i); if (g < 0.1 || g > 2 * PI - 0.1) return false;
  }
  const cr = (a, b, c) => (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  for (let i = 0; i < n; i++) for (let j = i + 2; j < n; j++) {
    if (i === 0 && j === n - 1) continue;
    const a = V(i), b = V(i + 1), c = V(j), d = V(j + 1);
    if (cr(a, b, c) * cr(a, b, d) < 0 && cr(c, d, a) * cr(c, d, b) < 0) return false;
  }
  return true;
}
