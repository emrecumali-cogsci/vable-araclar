/* V-Able Tek Denekli Araştırma Analiz Motoru (SCED engine)
 * Formüller SingleCaseES 0.7.4, scdhlm 0.7.4.9999 / lmeInfo ve scan 0.68.1 ile uyumlu olacak şekilde yazılmıştır.
 * MIT lisansı. */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.SCED = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /* ---------------- temel istatistik ---------------- */
  const sum = a => a.reduce((s, x) => s + x, 0);
  const mean = a => sum(a) / a.length;
  const variance = a => { if (a.length < 2) return NaN; const m = mean(a); return sum(a.map(x => (x - m) ** 2)) / (a.length - 1); };
  const sd = a => Math.sqrt(variance(a));
  const sorted = a => a.slice().sort((x, y) => x - y);
  function median(a) { const s = sorted(a), n = s.length; if (!n) return NaN; return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2; }

  // Wichura AS241 normal quantile
  function qnorm(p) {
    if (p <= 0) return -Infinity; if (p >= 1) return Infinity;
    const q = p - 0.5;
    let r, val;
    if (Math.abs(q) <= 0.425) {
      r = 0.180625 - q * q;
      val = q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r + 45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r + 133.14166789178437745) * r + 3.387132872796366608) /
        (((((((r * 5226.495278852545925 + 28729.085735721942674) * r + 39307.89580009271061) * r + 21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r + 42.313330701600911252) * r + 1);
      return val;
    }
    r = q < 0 ? p : 1 - p;
    r = Math.sqrt(-Math.log(r));
    if (r <= 5) {
      r -= 1.6;
      val = (((((((r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r + 0.24178072517745061177) * r + 1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) * r + 4.6303378461565452959) * r + 1.42343711074968357734) /
        (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r + 0.0151986665636164571966) * r + 0.14810397642748007459) * r + 0.68976733498510000455) * r + 1.6763848301838038494) * r + 2.05319162663775882187) * r + 1);
    } else {
      r -= 5;
      val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r + 0.0012426609473880784386) * r + 0.026532189526576123093) * r + 0.29656057182850489123) * r + 1.7848265399172913358) * r + 5.4637849111641143699) * r + 6.6579046435011037772) /
        (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r + 1.8463183175100546818e-5) * r + 7.868691311456132591e-4) * r + 0.0148753612908506148525) * r + 0.13692988092273580531) * r + 0.59983220655588793769) * r + 1);
    }
    return q < 0 ? -val : val;
  }
  // regularized incomplete beta (Numerical Recipes continued fraction) -> pt / qt
  function lgamma(x) {
    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, t = x + 5.5; t -= (x + 0.5) * Math.log(t); let s = 1.000000000190015;
    for (let j = 0; j < 6; j++) s += c[j] / ++y;
    return -t + Math.log(2.5066282746310005 * s / x);
  }
  function betacf(a, b, x) {
    const MAXIT = 300, EPS = 3e-16, FPMIN = 1e-300;
    let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d; let h = d;
    for (let m = 1; m <= MAXIT; m++) {
      const m2 = 2 * m; let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d;
      const del = d * c; h *= del; if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }
  function ibeta(x, a, b) {
    if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b;
  }
  function pt(t, df) { const x = df / (df + t * t); const p = 0.5 * ibeta(x, df / 2, 0.5); return t > 0 ? 1 - p : p; }
  function qt(p, df) {
    if (!isFinite(df) || df > 1e7) return qnorm(p);
    let lo = -1e3, hi = 1e3;
    for (let i = 0; i < 200; i++) { const mid = (lo + hi) / 2; if (pt(mid, df) < p) lo = mid; else hi = mid; if (hi - lo < 1e-12) break; }
    return (lo + hi) / 2;
  }
  function uniroot(f, a, b, tol = 1e-10) { // Brent
    let fa = f(a), fb = f(b);
    if (fa === 0) return a; if (fb === 0) return b;
    if (fa * fb > 0) return NaN;
    let c = a, fc = fa, d = b - a, e = d;
    for (let i = 0; i < 200; i++) {
      if (fb * fc > 0) { c = a; fc = fa; d = e = b - a; }
      if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
      const tol1 = 2 * 2.2e-16 * Math.abs(b) + tol / 2, xm = (c - b) / 2;
      if (Math.abs(xm) <= tol1 || fb === 0) return b;
      if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
        let s = fb / fa, p, q;
        if (a === c) { p = 2 * xm * s; q = 1 - s; }
        else { q = fa / fc; const r = fb / fc; p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1)); q = (q - 1) * (r - 1) * (s - 1); }
        if (p > 0) q = -q; p = Math.abs(p);
        if (2 * p < Math.min(3 * xm * q - Math.abs(tol1 * q), Math.abs(e * q))) { e = d; d = p / q; } else { d = xm; e = d; }
      } else { d = xm; e = d; }
      a = b; fa = fb;
      b += Math.abs(d) > tol1 ? d : (xm > 0 ? tol1 : -tol1);
      fb = f(b);
    }
    return b;
  }

  /* ---------------- örtüşme ölçüleri (SingleCaseES) ---------------- */
  const flip = (A, B, improvement) => improvement === "decrease" ? [A.map(x => -x), B.map(x => -x)] : [A, B];

  function NAP(A0, B0, improvement = "increase", confidence = 0.95) {
    const [A, B] = flip(A0, B0, improvement), m = A.length, n = B.length;
    const Q = A.map(a => B.map(b => (b > a ? 1 : 0) + (b === a ? 0.5 : 0)));
    const est = mean(Q.flat());
    const rowS = Q.map(r => sum(r.map(q => q - est)));
    const colS = B.map((_, j) => sum(A.map((_, i) => Q[i][j] - est)));
    const Q1 = sum(rowS.map(x => x * x)) / (m * n * n);
    const Q2 = sum(colS.map(x => x * x)) / (m * m * n);
    const tr = 0.5 / (m * n), Nt = Math.min(Math.max(est, tr), 1 - tr);
    const X = sum(Q.flat().map(q => (q - est) ** 2)) / (m * n);
    const V = (Nt * (1 - Nt) + n * Q1 + m * Q2 - 2 * X) / ((m - 1) * (n - 1));
    const h = (m + n) / 2 - 1, z = qnorm(1 - (1 - confidence) / 2);
    const f = x => m * n * (est - x) ** 2 * (2 - x) * (1 + x) - z * z * x * (1 - x) * (2 + h + (1 + 2 * h) * x * (1 - x));
    const lo = est > 0 ? uniroot(f, 0, est) : 0, hi = est < 1 ? uniroot(f, est, 1) : 1;
    return { ES: "NAP", Est: est, SE: Math.sqrt(V), CI_lower: lo, CI_upper: hi };
  }
  function Tau(A, B, improvement = "increase", confidence = 0.95) {
    const r = NAP(A, B, improvement, confidence);
    return { ES: "Tau", Est: 2 * r.Est - 1, SE: 2 * r.SE, CI_lower: 2 * r.CI_lower - 1, CI_upper: 2 * r.CI_upper - 1 };
  }
  function TauU(A0, B0, improvement = "increase") {
    const [A, B] = flip(A0, B0, improvement), m = A.length, n = B.length;
    let sP = 0; for (const b of B) for (const a of A) sP += (b > a) - (b < a);
    let sB = 0; for (let i = 0; i < m; i++) for (let j = i + 1; j < m; j++) sB += (A[j] > A[i]) - (A[j] < A[i]);
    return { ES: "Tau-U", Est: (sP - sB) / (m * n) };
  }
  function theilSen(y, x) {
    const s = []; for (let i = 0; i < y.length; i++) for (let j = i + 1; j < y.length; j++) s.push((y[j] - y[i]) / (x[j] - x[i]));
    const slope = median(s); return { slope, intercept: median(y.map((v, i) => v - x[i] * slope)) };
  }
  function TauBC(A, B, improvement = "increase", confidence = 0.95) {
    const m = A.length, n = B.length, sA = A.map((_, i) => i + 1), sB = B.map((_, i) => m + i + 1);
    const { slope, intercept } = theilSen(A, sA);
    const A2 = A.map((v, i) => v - intercept - slope * sA[i]), B2 = B.map((v, i) => v - intercept - slope * sB[i]);
    const r = Tau(A2, B2, improvement, confidence); r.ES = "Tau-BC"; r.slope = slope; r.intercept = intercept; return r;
  }
  function PND(A0, B0, improvement = "increase") { const [A, B] = flip(A0, B0, improvement), mx = Math.max(...A); return { ES: "PND", Est: mean(B.map(b => b > mx ? 1 : 0)) }; }
  function PEM(A0, B0, improvement = "increase") { const [A, B] = flip(A0, B0, improvement), md = median(A); return { ES: "PEM", Est: mean(B.map(b => (b > md ? 1 : 0) + (b === md ? 0.5 : 0))) }; }
  function PAND(A0, B0, improvement = "increase") {
    const [A, B] = flip(A0, B0, improvement), m = A.length, n = B.length;
    const X = [-Infinity, ...sorted(A)], Y = [...sorted(B), Infinity];
    let best = 0; for (let i = 1; i <= m + 1; i++) for (let j = 1; j <= n + 1; j++) if (X[i - 1] < Y[j - 1]) best = Math.max(best, i + n - j);
    return { ES: "PAND", Est: best / (m + n) };
  }
  function IRD(A, B, improvement = "increase") {
    const p = PAND(A, B, improvement).Est, m = A.length, n = B.length;
    return { ES: "IRD", Est: ((m + n) ** 2 * p - m * m - n * n) / (2 * m * n) };
  }

  /* ---------------- ortalama farkına dayalı ölçüler ---------------- */
  function truncConst(scale, intervals, obsLength) {
    if (scale === "count") return 1;
    if (scale === "rate") return obsLength > 0 ? obsLength : Infinity;
    if (scale === "proportion") return intervals > 0 ? intervals : Infinity;
    if (scale === "percentage") return intervals > 0 ? intervals / 100 : Infinity;
    return Infinity;
  }
  function LRR(A, B, kind, improvement, opts = {}) { // kind: "i" | "d"
    const scale = opts.scale || "count", D = opts.D_const ?? truncConst(scale, opts.intervals, opts.observation_length), conf = opts.confidence ?? 0.95;
    const n = [A.length, B.length]; let M = [mean(A), mean(B)], V = [variance(A), variance(B)];
    const other = kind === "i" ? "decrease" : "increase";
    if (improvement === other) { if (scale === "percentage") M = M.map(x => 100 - x); if (scale === "proportion") M = M.map(x => 1 - x); }
    if (D > 0) { M = M.map((x, i) => Math.max(1 / (2 * D * n[i]), x)); V = V.map((x, i) => Math.max(1 / (D * D * n[i] ** 3), x)); }
    const BC = M.map((x, i) => Math.log(x) + V[i] / (2 * n[i] * x * x));
    let est = BC[1] - BC[0];
    if (improvement === other && !["percentage", "proportion"].includes(scale)) est = -est;
    const SE = Math.sqrt(sum(V.map((v, i) => v / (n[i] * M[i] ** 2)))), z = qnorm(1 - (1 - conf) / 2);
    const pc = x => 100 * (Math.exp(x) - 1);
    return { ES: kind === "i" ? "LRRi" : "LRRd", Est: est, SE, CI_lower: est - z * SE, CI_upper: est + z * SE, pct: pc(est), pct_lower: pc(est - z * SE), pct_upper: pc(est + z * SE) };
  }
  const LRRi = (A, B, improvement = "increase", o) => LRR(A, B, "i", improvement, o);
  const LRRd = (A, B, improvement = "decrease", o) => LRR(A, B, "d", improvement, o);
  function SMD(A, B, improvement = "increase", std_dev = "baseline", confidence = 0.95) {
    const n = [A.length, B.length], M = [mean(A), mean(B)], V = [variance(A), variance(B)];
    let df, s2, SV1;
    if (std_dev === "baseline") { df = n[0] - 1; s2 = V[0]; SV1 = 1 / n[0] + V[1] / (n[1] * V[0]); }
    else { df = n[0] + n[1] - 2; s2 = ((n[0] - 1) * V[0] + (n[1] - 1) * V[1]) / df; SV1 = 1 / n[0] + 1 / n[1]; }
    const J = 1 - 3 / (4 * df - 1); let d = J * (M[1] - M[0]) / Math.sqrt(s2); if (improvement === "decrease") d = -d;
    const SE = J * Math.sqrt(SV1 + d * d / (2 * df)), z = qnorm(1 - (1 - confidence) / 2);
    return { ES: "SMD", Est: d, SE, CI_lower: d - z * SE, CI_upper: d + z * SE, SD: Math.sqrt(s2), std_dev };
  }

  /* ---------------- görsel analiz yardımcıları ---------------- */
  function olsSlope(y, x) { const mx = mean(x), my = mean(y), sxx = sum(x.map(v => (v - mx) ** 2)); return sxx > 0 ? { slope: sum(x.map((v, i) => (v - mx) * (y[i] - my))) / sxx, intercept: my - (sxx > 0 ? sum(x.map((v, i) => (v - mx) * (y[i] - my))) / sxx : 0) * mx } : { slope: NaN, intercept: my }; }
  // split-middle (White & Haring, 1980): veriler iki yarıya bölünür; tek sayıda ise orta nokta dışarıda bırakılır.
  function splitMiddle(y, x) {
    const n = y.length; if (n < 2) return { slope: NaN, intercept: NaN };
    const h = Math.floor(n / 2), y1 = y.slice(0, h), x1 = x.slice(0, h), y2 = y.slice(n - h), x2 = x.slice(n - h);
    const p1 = [median(x1), median(y1)], p2 = [median(x2), median(y2)];
    const slope = p2[0] === p1[0] ? NaN : (p2[1] - p1[1]) / (p2[0] - p1[0]);
    return { slope, intercept: p1[1] - slope * p1[0], p1, p2 };
  }
  function phaseSummary(y, x) {
    const md = median(y), ols = olsSlope(y, x), sm = splitMiddle(y, x);
    const within = md === 0 ? NaN : mean(y.map(v => Math.abs(v - md) <= 0.25 * Math.abs(md) ? 1 : 0));
    return { n: y.length, mean: mean(y), median: md, sd: sd(y), min: Math.min(...y), max: Math.max(...y), range: Math.max(...y) - Math.min(...y),
      ols_slope: ols.slope, ols_intercept: ols.intercept, sm_slope: sm.slope, sm_intercept: sm.intercept, stability: within, stable: within >= 0.8 };
  }
  function immediacy(A, B, k = 3) { const a = A.slice(-Math.min(k, A.length)), b = B.slice(0, Math.min(k, B.length)); return { lastA: mean(a), firstB: mean(b), diff: mean(b) - mean(a), k }; }

  /* ---------------- küçük matris araçları ---------------- */
  const zeros = (r, c) => Array.from({ length: r }, () => new Array(c).fill(0));
  function matmul(A, B) { const r = A.length, k = B.length, c = B[0].length, C = zeros(r, c); for (let i = 0; i < r; i++) for (let t = 0; t < k; t++) { const a = A[i][t]; if (a) for (let j = 0; j < c; j++) C[i][j] += a * B[t][j]; } return C; }
  const T = A => A[0].map((_, j) => A.map(r => r[j]));
  function chol(A) { const n = A.length, L = zeros(n, n); for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) { let s = A[i][j]; for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k]; if (i === j) { if (s <= 0) return null; L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j]; } return L; }
  function cholInv(A) { const L = chol(A); if (!L) return null; const n = A.length, Li = zeros(n, n); for (let i = 0; i < n; i++) { Li[i][i] = 1 / L[i][i]; for (let j = 0; j < i; j++) { let s = 0; for (let k = j; k < i; k++) s -= L[i][k] * Li[k][j]; Li[i][j] = s / L[i][i]; } } const inv = matmul(T(Li), Li); let logdet = 0; for (let i = 0; i < n; i++) logdet += 2 * Math.log(L[i][i]); return { inv, logdet }; }
  const trProd = (A, B) => { let s = 0; for (let i = 0; i < A.length; i++) for (let j = 0; j < A.length; j++) s += A[i][j] * B[j][i]; return s; };
  function nelderMead(f, x0, step = 0.5, iters = 2000, tol = 1e-12) {
    const n = x0.length; let P = [x0.slice()]; for (let i = 0; i < n; i++) { const p = x0.slice(); p[i] += step; P.push(p); }
    let F = P.map(f);
    for (let it = 0; it < iters; it++) {
      const idx = F.map((v, i) => i).sort((a, b) => F[a] - F[b]); P = idx.map(i => P[i]); F = idx.map(i => F[i]);
      if (Math.abs(F[n] - F[0]) < tol * (1 + Math.abs(F[0]))) break;
      const c = new Array(n).fill(0); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) c[j] += P[i][j] / n;
      const xr = c.map((v, j) => v + (v - P[n][j])), fr = f(xr);
      if (fr < F[0]) { const xe = c.map((v, j) => v + 2 * (v - P[n][j])), fe = f(xe); if (fe < fr) { P[n] = xe; F[n] = fe; } else { P[n] = xr; F[n] = fr; } }
      else if (fr < F[n - 1]) { P[n] = xr; F[n] = fr; }
      else { const xc = c.map((v, j) => v + 0.5 * (P[n][j] - v)), fc = f(xc); if (fc < F[n]) { P[n] = xc; F[n] = fc; } else { for (let i = 1; i <= n; i++) { P[i] = P[i].map((v, j) => P[0][j] + 0.5 * (v - P[0][j])); F[i] = f(P[i]); } } }
    }
    const b = F.indexOf(Math.min(...F)); return { x: P[b], f: F[b] };
  }

  /* ---------------- BC-SMD: çoklu başlama düzeyi, temel model ----------------
   * y_ij = b0 + b1*Tedavi_ij + u_j + e_ij, u_j ~ N(0, tau^2), e_ij AR(1) (phi, sigma^2), REML.
   * scdhlm::calc_BCSMD(design = "MBP", FE_base = 0, RE_base = 0, FE_trt = 0, corStruct = "AR1") ile eşdeğerdir. */
  function BCSMD_MB(cases, confidence = 0.95) {
    // cases: [{ y:[], trt:[0/1], time:[] }]
    const blocks = cases.map(c => ({ y: c.y, X: c.y.map((_, i) => [1, c.trt[i]]), t: c.time }));
    const N = sum(blocks.map(b => b.y.length)), p = 2;
    const Rmat = (t, phi) => t.map(a => t.map(b => phi ** Math.abs(a - b)));
    function parts(lam, phi) { // W = lam*J + R
      let XtWX = zeros(p, p), XtWy = [0, 0], logdetW = 0; const Winv = [];
      for (const b of blocks) {
        const R = Rmat(b.t, phi), W = R.map((r, i) => r.map((v, j) => v + lam));
        const ci = cholInv(W); if (!ci) return null; logdetW += ci.logdet; Winv.push(ci.inv);
        const WX = matmul(ci.inv, b.X), XtWXb = matmul(T(b.X), WX);
        for (let i = 0; i < p; i++) { for (let j = 0; j < p; j++) XtWX[i][j] += XtWXb[i][j]; XtWy[i] += sum(b.y.map((y, k) => WX[k][i] * y)); }
      }
      const ci = cholInv(XtWX); if (!ci) return null;
      const beta = [0, 1].map(i => ci.inv[i][0] * XtWy[0] + ci.inv[i][1] * XtWy[1]);
      let rss = 0; blocks.forEach((b, k) => { const r = b.y.map((y, i) => y - beta[0] * b.X[i][0] - beta[1] * b.X[i][1]); for (let i = 0; i < r.length; i++) for (let j = 0; j < r.length; j++) rss += r[i] * Winv[k][i][j] * r[j]; });
      return { beta, rss, logdetW, logdetXtWX: ci.logdet, XtWXinv: ci.inv, Winv };
    }
    const negll = th => { const lam = Math.exp(th[0]), phi = Math.tanh(th[1]); const q = parts(lam, phi); if (!q) return 1e100; return 0.5 * ((N - p) * Math.log(q.rss / (N - p)) + q.logdetW + q.logdetXtWX); };
    let best = null;
    for (const s of [[Math.log(0.5), 0.2], [Math.log(0.05), 0], [Math.log(2), 0.5], [Math.log(1e-3), -0.2], [Math.log(5), 0.8]]) {
      const r = nelderMead(negll, s, 0.6, 4000, 1e-14); if (!best || r.f < best.f) best = r;
    }
    let r2 = nelderMead(negll, best.x, 0.1, 4000, 1e-15); if (r2.f < best.f) best = r2;
    const lam = Math.exp(best.x[0]), phi = Math.tanh(best.x[1]), q = parts(lam, phi), sigma2 = q.rss / (N - p), tau2 = lam * sigma2;
    const beta = q.beta, vcovBeta = q.XtWXinv.map(r => r.map(v => v * sigma2));
    // beklenen REML bilgi matrisi: I_ij = 0.5 tr(P dV_i P dV_j), theta = (tau^2, phi, sigma^2)
    const Vinv = q.Winv.map(W => W.map(r => r.map(v => v / sigma2)));
    const Xall = [], dV = [[], [], []]; let offset = 0;
    const XtVinvX_inv = q.XtWXinv.map(r => r.map(v => v * sigma2));
    // P = V^-1 - V^-1 X (X'V^-1X)^-1 X'V^-1 (tam boyutta; N küçük)
    const Vfull = zeros(N, N), Xf = [];
    blocks.forEach((b, k) => { const n = b.y.length; for (let i = 0; i < n; i++) { Xf.push(b.X[i]); for (let j = 0; j < n; j++) Vfull[offset + i][offset + j] = Vinv[k][i][j]; } offset += n; });
    const VX = matmul(Vfull, Xf), Pm = (() => { const A = matmul(matmul(VX, XtVinvX_inv), T(VX)); return Vfull.map((r, i) => r.map((v, j) => v - A[i][j])); })();
    const D = [zeros(N, N), zeros(N, N), zeros(N, N)]; offset = 0;
    blocks.forEach(b => { const n = b.t.length; for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { const d = Math.abs(b.t[i] - b.t[j]);
      D[0][offset + i][offset + j] = 1; D[1][offset + i][offset + j] = d === 0 ? 0 : sigma2 * d * phi ** (d - 1); D[2][offset + i][offset + j] = phi ** d; } offset += n; });
    const PD = D.map(d => matmul(Pm, d));
    const I = zeros(3, 3); for (let i = 0; i < 3; i++) for (let j = 0; j <= i; j++) I[i][j] = I[j][i] = 0.5 * trProd(PD[i], PD[j]);
    const Iinv = cholInv(I);
    const rc = [1, 0, 1], r_theta = tau2 + sigma2, delta = beta[1] / Math.sqrt(r_theta);
    let varR = 0; for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) varR += rc[i] * rc[j] * (Iinv ? Iinv.inv[i][j] : NaN);
    const nu = 2 * r_theta ** 2 / varR, J = 1 - 3 / (4 * nu - 1), g = J * delta, kappa2 = vcovBeta[1][1] / r_theta;
    const nut = Math.max(nu, 2.001), Jt = 1 - 3 / (4 * nut - 1);
    const Vg = J * J * (nut * kappa2 / (nut - 2) + g * g * (nut / (nut - 2) - 1 / (Jt * Jt))), SEg = Math.sqrt(Vg);
    const tq = qt(1 - (1 - confidence) / 2, nu);
    return { ES: "BC-SMD", Est: g, SE: SEg, CI_lower: g - tq * SEg, CI_upper: g + tq * SEg, delta, nu, J, kappa: Math.sqrt(kappa2),
      beta, tau2, sigma2, phi, rho: tau2 / r_theta, logLik: -best.f - 0.5 * (N - p) * (1 + Math.log(2 * Math.PI)) };
  }

  /* ---------------- veri hazırlama ---------------- */
  // kayıtlar: {case, phase, session, outcome, condition?, criterion?}
  function phaseRuns(rows) { // ardışık aynı etiketler = evre
    const runs = []; let cur = null;
    for (const r of rows) { if (!cur || cur.label !== r.phase) { cur = { label: r.phase, rows: [] }; runs.push(cur); } cur.rows.push(r); }
    const baseType = runs.length ? phaseType(runs[0].label) : "A";
    runs.forEach((r, i) => { r.type = phaseType(r.label); r.isBaseline = r.type === baseType; r.index = i; });
    return runs;
  }
  function phaseType(label) { const s = String(label).trim().toUpperCase(); const m = s.match(/^[A-ZÇĞİÖŞÜ]+/); return m ? m[0] : s; }
  function comparisons(runs) { // her başlama evresi ardından gelen müdahale evresiyle eşleşir (A1-B1, A2-B2 ...)
    const out = []; for (let i = 0; i < runs.length - 1; i++) if (runs[i].isBaseline && !runs[i + 1].isBaseline) out.push([runs[i], runs[i + 1]]); return out;
  }

  const ALL = { NAP, Tau, TauU, TauBC, PND, PEM, PAND, IRD, LRRi, LRRd, SMD };
  return { mean, median, variance, sd, qnorm, qt, pt, uniroot, theilSen, ...ALL, olsSlope, splitMiddle, phaseSummary, immediacy, BCSMD_MB, phaseRuns, phaseType, comparisons };
});
