// JS motoru ile R referanslarını karşılaştırır
const S = require("./engine.js");
const ref = require("./ref.json");
const out = { rows: [], maxdiff: {} };
const tolES = 1e-9;
function cmp(label, field, js, r, tol) {
  // NAP/Tau GA kökleri: R uniroot varsayılan toleransı ~1.2e-4 (bkz. ci_tight.R)
  if (/^CI_/.test(field) && ["NAP", "Tau", "TauBC"].includes(label)) tol = 2e-4;
  if (r === null || r === undefined || (typeof r === "number" && !isFinite(r)) && !isFinite(js)) return;
  const d = Math.abs(js - r);
  const key = label + "." + field;
  out.maxdiff[key] = Math.max(out.maxdiff[key] || 0, isNaN(d) ? Infinity : d);
  if (!(d <= tol)) out.rows.push({ key, js, r, d });
}
const first = x => Array.isArray(x) ? x[0] : x;
const get = (o, f) => { o = first(o); return o ? first(o[f]) : undefined; };
for (const p of ref.pairs) {
  const A = [].concat(p.A), B = [].concat(p.B), imp = p.improvement;
  const run = (name, js, fields) => { const r = p[name]; if (!r) return; for (const f of fields) cmp(name, f, js[f], get(r, f), tolES); };
  run("NAP", S.NAP(A, B, imp), ["Est", "SE", "CI_lower", "CI_upper"]);
  run("Tau", S.Tau(A, B, imp), ["Est", "SE", "CI_lower", "CI_upper"]);
  run("TauU", S.TauU(A, B, imp), ["Est"]);
  run("TauBC", S.TauBC(A, B, imp), ["Est", "SE", "CI_lower", "CI_upper"]);
  run("PND", S.PND(A, B, imp), ["Est"]);
  run("PEM", S.PEM(A, B, imp), ["Est"]);
  run("PAND", S.PAND(A, B, imp), ["Est"]);
  run("IRD", S.IRD(A, B, imp), ["Est"]);
  run("SMDb", S.SMD(A, B, imp, "baseline"), ["Est", "SE", "CI_lower", "CI_upper"]);
  run("SMDp", S.SMD(A, B, imp, "pool"), ["Est", "SE", "CI_lower", "CI_upper"]);
  for (const k of ["LRRi", "LRRd"]) {
    const r = p[k]; if (!r) continue;
    const js = S[k](A, B, imp, { scale: "count" });
    cmp(k, "Est", js.Est, r[0].Est, tolES); cmp(k, "SE", js.SE, r[0].SE, tolES);
    cmp(k, "CI_lower", js.CI_lower, r[0].CI_lower, tolES); cmp(k, "pct", js.pct, r[1].Est, 1e-7); cmp(k, "pct_lower", js.pct_lower, r[1].CI_lower, 1e-7);
  }
}
for (const b of ref.bcsmd) {
  if (b.ref.error) { out.rows.push({ key: "BCSMD." + b.name, err: b.ref.error }); continue; }
  const byCase = {};
  for (const row of b.data) { const c = row.case; (byCase[c] = byCase[c] || { y: [], trt: [], time: [] }); byCase[c].y.push(row.outcome); byCase[c].trt.push(row.trt); byCase[c].time.push(row.session); }
  const t0 = Date.now(); const js = S.BCSMD_MB(Object.values(byCase)); const ms = Date.now() - t0;
  const r = b.ref, rel = (a, c) => Math.abs(a - c) / Math.max(1, Math.abs(c));
  const line = { name: b.name, ms, g: [js.Est, r.g], SE: [js.SE, r.SE], nu: [js.nu, r.nu], phi: [js.phi, r.phi], tau2: [js.tau2, r.tau2], sigma2: [js.sigma2, r.sigma2], logLik: [js.logLik, r.logLik], CI: [[js.CI_lower, js.CI_upper], [r.CI_L, r.CI_U]] };
  line.rel = { g: rel(js.Est, r.g), SE: rel(js.SE, r.SE), nu: rel(js.nu, r.nu), CI_L: rel(js.CI_lower, r.CI_L), CI_U: rel(js.CI_upper, r.CI_U), logLik: js.logLik - r.logLik };
  (out.bc = out.bc || []).push(line);
}
console.log("Örtüşme/parametrik — en büyük mutlak fark (tolerans 1e-9):");
for (const [k, v] of Object.entries(out.maxdiff)) console.log("  " + k.padEnd(18), v.toExponential(2));
console.log("Tolerans dışı:", out.rows.length); out.rows.slice(0, 15).forEach(r => console.log("  ", JSON.stringify(r)));
console.log("\nBC-SMD (göreli fark; logLik mutlak fark):");
for (const l of out.bc || []) console.log("  " + l.name.padEnd(17), Object.entries(l.rel).map(([k, v]) => k + "=" + v.toExponential(1)).join(" "), " g=" + l.g.map(x => x.toFixed(5)).join("/"), l.ms + "ms");
require("fs").writeFileSync(__dirname + "/compare_out.json", JSON.stringify(out, null, 1));
