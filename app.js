/* Fluidized Beds — UI logic.
   Physics comes from build/release.wasm (AssemblyScript). If the WASM
   module fails to load, an identical plain-JS implementation is used as a
   fallback so the site still works. */

// ---------- Presets ----------
const GASES = {
  air20:  { name: "Air, 20 °C",      rho: 1.20, mu: 1.81e-5, tempC: 20 },
  air500: { name: "Hot air, 500 °C", rho: 0.46, mu: 3.70e-5, tempC: 500 },
  n2:     { name: "Nitrogen, 20 °C", rho: 1.16, mu: 1.76e-5, tempC: 20 },
  steam:  { name: "Steam, 200 °C",   rho: 0.59, mu: 1.60e-5, tempC: 200 },
};
const PARTICLES = {
  sand:    { name: "Sand",         rhos: 2600, dpUm: 500 },
  glass:   { name: "Glass beads",  rhos: 2500, dpUm: 500 },
  fcc:     { name: "FCC catalyst", rhos: 1500, dpUm: 75 },
  alumina: { name: "Alumina",      rhos: 3900, dpUm: 500 },
  coal:    { name: "Coal",         rhos: 1400, dpUm: 1000 },
};

// Air properties at 1 atm as a function of temperature (ideal gas + Sutherland-ish fit).
function airProps(tempC) {
  const T = tempC + 273.15;
  const rho = 1.20 * (293.15 / T);
  const mu = 1.81e-5 * Math.pow(T / 293.15, 0.7);
  return { rho, mu };
}

// ---------- Shared state (defaults: air 20 C, sand 500 um, H=0.5 m, eps=0.4) ----------
const S = {
  dp: 500e-6, rhos: 2600, rhog: 1.20, mu: 1.81e-5,
  eps: 0.40, H: 0.50, tempC: 20, gas: "air20",
};

// ---------- Physics: WASM first, JS fallback ----------
const JSPhys = {
  ergunDeltaPPerL: (u, mu, rho, eps, dp) => {
    const e3 = eps ** 3;
    return 150 * mu * (1 - eps) ** 2 * u / (e3 * dp * dp)
         + 1.75 * (1 - eps) * rho * u * u / (e3 * dp);
  },
  ergunViscousPerL: (u, mu, eps, dp) => 150 * mu * (1 - eps) ** 2 * u / (eps ** 3 * dp * dp),
  ergunInertialPerL: (u, rho, eps, dp) => 1.75 * (1 - eps) * rho * u * u / (eps ** 3 * dp),
  bedWeightDeltaP: (eps, rhos, rhog, H) => (1 - eps) * (rhos - rhog) * 9.81 * H,
  umfErgun: (mu, rho, eps, dp, rhos, rhog) => {
    const e3 = eps ** 3;
    const A = 150 * mu * (1 - eps) ** 2 / (e3 * dp * dp);
    const B = 1.75 * (1 - eps) * rho / (e3 * dp);
    const C = (1 - eps) * (rhos - rhog) * 9.81;
    return (-A + Math.sqrt(A * A + 4 * B * C)) / (2 * B);
  },
  umfSmallParticle: (mu, eps, dp, rhos, rhog) =>
    dp * dp * (rhos - rhog) * 9.81 * eps ** 3 / (150 * mu * (1 - eps)),
  archimedes: (dp, rhos, rhog, mu) => dp ** 3 * rhog * (rhos - rhog) * 9.81 / (mu * mu),
  umfWenYu: (dp, rhos, rhog, mu) => {
    const Ar = dp ** 3 * rhog * (rhos - rhog) * 9.81 / (mu * mu);
    return (Math.sqrt(33.7 * 33.7 + 0.0408 * Ar) - 33.7) * mu / (dp * rhog);
  },
  dragCoefficient: (Re) => Re < 1000 ? (24 / Re) * (1 + 0.15 * Math.pow(Re, 0.687)) : 0.44,
  terminalVelocity: (dp, rhos, rhog, mu) => {
    let ut = dp * dp * (rhos - rhog) * 9.81 / (18 * mu);
    for (let i = 0; i < 200; i++) {
      const Re = rhog * ut * dp / mu;
      const Cd = Re < 1000 ? (24 / Re) * (1 + 0.15 * Math.pow(Re, 0.687)) : 0.44;
      const next = Math.sqrt(4 * dp * (rhos - rhog) * 9.81 / (3 * rhog * Cd));
      if (Math.abs(next - ut) < 1e-9 * next) { ut = next; break; }
      ut = next;
    }
    return ut;
  },
  reynoldsTerminal: (ut, dp, rhog, mu) => rhog * ut * dp / mu,
  richardsonZakiN: (ReT) => {
    if (ReT < 0.2) return 4.65;
    if (ReT < 1) return 4.4 * Math.pow(ReT, -0.03);
    if (ReT < 500) return 4.4 * Math.pow(ReT, -0.1);
    return 2.4;
  },
  voidageAtVelocity: (u, ut, ReT) => {
    const n = JSPhys.richardsonZakiN(ReT);
    const r = u / ut;
    if (r <= 0) return 0.4;
    if (r >= 1) return 1;
    return Math.pow(r, 1 / n);
  },
  bedHeight: (Hmf, epsMf, eps) => Hmf * (1 - epsMf) / (1 - eps),
  geldartGroup: (dp, dRho) => {
    const dpUm = dp * 1e6;
    if (dpUm < 30) return 0;
    const ref = Math.max(dRho, 1);
    const dpAB = 100 * Math.sqrt(1000 / ref);
    const dpBD = 1000 * Math.sqrt(1000 / ref);
    if (dpUm < dpAB) return 1;
    if (dpUm < dpBD) return 2;
    return 3;
  },
};

let P = JSPhys;       // active physics implementation
let wasmOn = false;

async function loadWasm() {
  try {
    const res = await fetch("build/release.wasm");
    if (!res.ok) throw new Error("HTTP " + res.status);
    const { instance } = await WebAssembly.instantiateStreaming(res, {});
    P = instance.exports;
    wasmOn = true;
    console.log("Fluidized Beds: WebAssembly physics loaded.");
  } catch (e) {
    console.warn("Fluidized Beds: WASM failed to load, using JS fallback.", e);
    P = JSPhys;
  }
  initAll();
}

// ---------- Small helpers ----------
const $ = (id) => document.getElementById(id);
const fmt = (x, d = 2) => (x == null || !isFinite(x)) ? "—" : x.toFixed(d);
const fmtSci = (x) => (x == null || !isFinite(x)) ? "—" : x.toExponential(2);
function regimeAt(u, umf, ut) {
  if (u < umf) return "Fixed bed";
  if (u < 3 * umf) return "Bubbling";
  if (u < ut) return "Turbulent";
  return "Transport (entrainment)";
}

// ---------- Tabs ----------
function initTabs() {
  const btns = document.querySelectorAll("nav.tabs button");
  btns.forEach((b) => b.addEventListener("click", () => showTab(b.dataset.tab)));
  const hash = location.hash.replace("#", "");
  if (hash && document.getElementById(hash)) showTab(hash);
}
function showTab(id) {
  document.querySelectorAll("nav.tabs button").forEach((b) => b.classList.toggle("active", b.dataset.tab === id));
  document.querySelectorAll(".tab-panel").forEach((p) => p.classList.toggle("active", p.id === id));
  history.replaceState(null, "", "#" + id);
  const panel = document.getElementById(id);
  if (panel) {
    // Redraw canvases in the newly shown tab (they need layout to be visible).
    requestAnimationFrame(() => { drawAll(); });
  }
}

// ---------- KaTeX ----------
function renderMath() {
  if (window.renderMathInElement) {
    renderMathInElement(document.body, {
      delimiters: [
        { left: "\\[", right: "\\]", display: true },
        { left: "\\(", right: "\\)", display: false },
      ],
      throwOnError: false,
    });
  }
}

// ---------- Canvas helper ----------
function setupCanvas(cv) {
  const dpr = window.devicePixelRatio || 1;
  const w = cv.clientWidth || cv.width;
  const h = (cv.height / cv.width) * w;
  cv.width = Math.round(w * dpr);
  cv.height = Math.round(h * dpr);
  const ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, w, h };
}
function drawAxes(ctx, x0, y0, x1, y1, xLabel, yLabel, xTicks, yTicks) {
  ctx.strokeStyle = "#555";
  ctx.fillStyle = "#a8a8ad";
  ctx.font = "12px sans-serif";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(x0, y0); ctx.lineTo(x1, y0); ctx.moveTo(x0, y0); ctx.lineTo(x0, y1);
  ctx.stroke();
  ctx.textAlign = "center";
  xTicks.forEach(([x, label]) => { ctx.fillText(label, x, y0 + 16); });
  ctx.textAlign = "right";
  yTicks.forEach(([y, label]) => { ctx.fillText(label, x0 - 6, y + 4); });
  ctx.fillStyle = "#e9e9e9";
  ctx.textAlign = "center";
  ctx.fillText(xLabel, (x0 + x1) / 2, y0 + 34);
  ctx.save();
  ctx.translate(14, (y0 + y1) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();
}

// ---------- Tab 3: Pressure-drop vs velocity graph ----------
function currentUmf() {
  return {
    ergun: P.umfErgun(S.mu, S.rhog, S.eps, S.dp, S.rhos, S.rhog),
    wenYu: P.umfWenYu(S.dp, S.rhos, S.rhog, S.mu),
    ut: P.terminalVelocity(S.dp, S.rhos, S.rhog, S.mu),
    dpmf: P.bedWeightDeltaP(S.eps, S.rhos, S.rhog, S.H),
  };
}

function drawDpGraph() {
  const cv = $("dpCanvas");
  if (!cv || !cv.closest(".tab-panel").classList.contains("active")) return;
  const { ctx, w, h } = setupCanvas(cv);
  const { ergun: umf, ut, dpmf } = currentUmf();
  const uMax = Math.min(Math.max(2.5 * umf, 0.5), 5);
  const yMax = dpmf * 1.25;
  const L = 70, R = 18, T = 16, B = 46;
  const X = (u) => L + (u / uMax) * (w - L - R);
  const Y = (dp) => h - B - (dp / yMax) * (h - T - B);

  // Regime bands.
  const bands = [
    [0, umf, "Fixed bed", "rgba(77,163,255,0.07)"],
    [umf, Math.min(3 * umf, ut), "Bubbling", "rgba(245,179,1,0.08)"],
    [Math.min(3 * umf, ut), ut, "Turbulent", "rgba(224,1,34,0.07)"],
    [ut, uMax, "Transport", "rgba(224,1,34,0.13)"],
  ];
  bands.forEach(([a, b, label, col]) => {
    if (b <= a || a >= uMax) return;
    ctx.fillStyle = col;
    ctx.fillRect(X(a), T, X(Math.min(b, uMax)) - X(a), h - B - T);
    ctx.fillStyle = "#a8a8ad";
    ctx.font = "11px sans-serif";
    ctx.textAlign = "center";
    const cx = (X(a) + X(Math.min(b, uMax))) / 2;
    if (X(Math.min(b, uMax)) - X(a) > 44) ctx.fillText(label, cx, T + 14);
  });

  // Ergun curve up to u_mf, then plateau.
  ctx.strokeStyle = "#e00122";
  ctx.lineWidth = 2.5;
  ctx.beginPath();
  const N = 200;
  for (let i = 0; i <= N; i++) {
    const u = (i / N) * uMax;
    const dpv = u <= umf ? P.ergunDeltaPPerL(u, S.mu, S.rhog, S.eps, S.dp) * S.H : dpmf;
    if (i === 0) ctx.moveTo(X(u), Y(dpv)); else ctx.lineTo(X(u), Y(dpv));
  }
  ctx.stroke();

  // Bed-weight plateau line (dashed).
  ctx.strokeStyle = "#f5b301";
  ctx.lineWidth = 1.5;
  ctx.setLineDash([6, 4]);
  ctx.beginPath();
  ctx.moveTo(X(0), Y(dpmf));
  ctx.lineTo(X(uMax), Y(dpmf));
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = "#f5b301";
  ctx.font = "12px sans-serif";
  ctx.textAlign = "left";
  ctx.fillText("ΔP = bed weight", X(uMax) - 132, Y(dpmf) - 8);

  // u_mf marker.
  ctx.strokeStyle = "#fff";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(X(umf), Y(0));
  ctx.lineTo(X(umf), Y(dpmf));
  ctx.stroke();
  ctx.fillStyle = "#e00122";
  ctx.beginPath();
  ctx.arc(X(umf), Y(dpmf), 6, 0, 2 * Math.PI);
  ctx.fill();
  ctx.fillStyle = "#fff";
  ctx.font = "bold 12px sans-serif";
  ctx.textAlign = "left";
  ctx.fillText("u_mf", X(umf) + 8, Y(dpmf) + 4);

  // ut marker.
  if (ut < uMax) {
    ctx.strokeStyle = "#4da3ff";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(X(ut), Y(0));
    ctx.lineTo(X(ut), Y(dpmf));
    ctx.stroke();
    ctx.fillStyle = "#4da3ff";
    ctx.font = "bold 12px sans-serif";
    ctx.fillText("u_t", X(ut) + 6, Y(dpmf) + 4);
  }

  const xTicks = [];
  for (let i = 0; i <= 5; i++) {
    const u = (i / 5) * uMax;
    xTicks.push([X(u), u.toFixed(2)]);
  }
  const yTicks = [];
  for (let i = 0; i <= 4; i++) {
    const dpv = (i / 4) * yMax;
    yTicks.push([Y(dpv), dpv >= 1000 ? (dpv / 1000).toFixed(1) + "k" : dpv.toFixed(0)]);
  }
  drawAxes(ctx, L, h - B, w - R, T, "Superficial gas velocity  u  (m/s)", "Pressure drop  ΔP  (Pa)", xTicks, yTicks);

  // Readouts.
  $("roUmfErgun").textContent = fmt(umf, 3) + " m/s";
  $("roUmfWY").textContent = fmt(currentUmf().wenYu, 3) + " m/s";
  $("roDPmf").textContent = fmt(dpmf, 0) + " Pa";
  $("roUt").textContent = fmt(ut, 2) + " m/s";
}

function bindDpControls() {
  const upd = () => {
    S.dp = parseFloat($("dpSlider").value) * 1e-6;
    S.rhos = parseFloat($("rhosSlider").value);
    S.rhog = parseFloat($("rhogSlider").value);
    S.mu = parseFloat($("muSlider").value) * 1e-5;
    S.eps = parseFloat($("epsSlider").value);
    S.H = parseFloat($("hSlider").value);
    $("dpVal").textContent = $("dpSlider").value + " μm";
    $("rhosVal").textContent = $("rhosSlider").value + " kg/m³";
    $("rhogVal").textContent = fmt(S.rhog, 2) + " kg/m³";
    $("muVal").innerHTML = fmt(parseFloat($("muSlider").value), 2) + " × 10<sup>−5</sup> Pa·s";
    $("epsVal").textContent = fmt(S.eps, 2);
    $("hVal").textContent = fmt(S.H, 2) + " m";
    drawDpGraph();
    drawBedFrame();
    updateCalc();
  };
  ["dpSlider", "rhosSlider", "rhogSlider", "muSlider", "epsSlider", "hSlider"].forEach((id) =>
    $(id).addEventListener("input", upd));

  $("tempSlider").addEventListener("input", () => {
    const t = parseFloat($("tempSlider").value);
    S.tempC = t;
    const ap = airProps(t);
    S.rhog = ap.rho;
    S.mu = ap.mu;
    S.gas = "air";
    $("tempVal").textContent = t + " °C";
    $("rhogSlider").value = S.rhog.toFixed(2);
    $("muSlider").value = (S.mu * 1e5).toFixed(2);
    document.querySelectorAll("#gasButtons .chip").forEach((c) => c.classList.remove("active"));
    upd();
  });

  document.querySelectorAll("#gasButtons .chip").forEach((btn) => {
    btn.addEventListener("click", () => {
      const g = GASES[btn.dataset.gas];
      S.gas = btn.dataset.gas;
      S.rhog = g.rho;
      S.mu = g.mu;
      S.tempC = g.tempC;
      $("rhogSlider").value = g.rho.toFixed(2);
      $("muSlider").value = (g.mu * 1e5).toFixed(2);
      $("tempSlider").value = g.tempC;
      $("tempVal").textContent = g.tempC + " °C";
      document.querySelectorAll("#gasButtons .chip").forEach((c) => c.classList.toggle("active", c === btn));
      upd();
    });
  });
  upd();
}

// ---------- Bed animation ----------
let animT = 0;
let particles = [];
function initParticles() {
  particles = [];
  const cols = 9, rows = 13;
  for (let r = 0; r < rows; r++) {
    for (let c = 0; c < cols; c++) {
      particles.push({
        x: 0.12 + c * 0.095 + (r % 2 ? 0.045 : 0),
        y: 0.86 - r * 0.058,
        ph: Math.random() * Math.PI * 2,
        sp: 0.5 + Math.random(),
      });
    }
  }
}
function drawBedFrame() {
  const cv = $("bedCanvas");
  if (!cv || !cv.closest(".tab-panel").classList.contains("active")) return;
  const { ctx, w, h } = setupCanvas(cv);
  const u = parseFloat($("animU").value);
  const { ergun: umf, ut } = currentUmf();
  const ratio = umf > 0 ? u / umf : 0;

  // Column.
  const cx0 = w * 0.18, cx1 = w * 0.82, top = 14, bot = h - 34;
  ctx.strokeStyle = "#e9e9e9";
  ctx.lineWidth = 2;
  ctx.strokeRect(cx0, top, cx1 - cx0, bot - top);
  // Distributor.
  ctx.strokeStyle = "#e00122";
  ctx.lineWidth = 4;
  ctx.beginPath();
  ctx.moveTo(cx0, bot - 26);
  ctx.lineTo(cx1, bot - 26);
  ctx.stroke();
  // Gas in arrow.
  ctx.strokeStyle = "#4da3ff";
  ctx.lineWidth = 2.5;
  ctx.beginPath();
  ctx.moveTo(w / 2, h - 4);
  ctx.lineTo(w / 2, bot - 18);
  ctx.stroke();

  // Particles.
  const fluid = u >= umf;
  const jig = fluid ? Math.min(1, (u - umf) / Math.max(umf, 0.05)) * 7 + 1.5 : 0;
  const spread = fluid ? Math.min(0.10, 0.02 + 0.05 * ratio) : 0;
  particles.forEach((p) => {
    const jx = fluid ? Math.sin(animT * 2 * p.sp + p.ph) * jig : 0;
    const jy = fluid ? Math.cos(animT * 1.7 * p.sp + p.ph) * jig * 0.7 : 0;
    const px = cx0 + (p.x + spread * (p.x - 0.5) * 2) * (cx1 - cx0) + jx;
    const py = top + (p.y - spread * 0.6) * (bot - 26 - top) + jy;
    ctx.fillStyle = "#e00122";
    ctx.beginPath();
    ctx.arc(px, py, 7, 0, 2 * Math.PI);
    ctx.fill();
  });

  // Bubbles when well above u_mf.
  if (u > 1.3 * umf && u < ut) {
    for (let i = 0; i < 3; i++) {
      const by = bot - 40 - ((animT * 40 * (1 + i * 0.3) + i * 90) % (bot - 70));
      const bx = cx0 + (cx1 - cx0) * (0.3 + 0.2 * i);
      ctx.strokeStyle = "#f5b301";
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(bx, by, 12 + i * 3, 0, 2 * Math.PI);
      ctx.stroke();
    }
  }
  // Entrainment: particles escaping above u_t.
  if (u >= ut) {
    ctx.fillStyle = "#e00122";
    for (let i = 0; i < 5; i++) {
      const ey = top - ((animT * 60 + i * 30) % 40) + 20;
      ctx.beginPath();
      ctx.arc(cx0 + (cx1 - cx0) * (0.2 + 0.15 * i), Math.max(4, ey), 5, 0, 2 * Math.PI);
      ctx.fill();
    }
  }

  $("animUVal").textContent = fmt(u, 2) + " m/s";
  $("animRegime").textContent = regimeAt(u, umf, ut);
  $("animRatio").textContent = fmt(ratio, 2);
}
function bedLoop() {
  animT += 0.05;
  drawBedFrame();
  requestAnimationFrame(bedLoop);
}

// ---------- u_mf calculator ----------
function updateCalc() {
  const dp = parseFloat($("calcDp").value) * 1e-6;
  const rhos = parseFloat($("calcRhos").value);
  const eps = parseFloat($("calcEps").value);
  const u = parseFloat($("calcU").value);
  const gas = GASES[$("calcGas").value];
  if (!isFinite(dp) || !isFinite(rhos) || !isFinite(eps)) return;
  const umfE = P.umfErgun(gas.mu, gas.rho, eps, dp, rhos, gas.rho);
  const umfW = P.umfWenYu(dp, rhos, gas.rho, gas.mu);
  const umfS = P.umfSmallParticle(gas.mu, eps, dp, rhos, gas.rho);
  const ut = P.terminalVelocity(dp, rhos, gas.rho, gas.mu);
  $("calcUmfErgun").textContent = fmt(umfE, 3) + " m/s";
  $("calcUmfWY").textContent = fmt(umfW, 3) + " m/s";
  $("calcUmfSmall").textContent = fmt(umfS, 3) + " m/s";
  $("calcUt").textContent = fmt(ut, 2) + " m/s";
  $("calcRegime").textContent = isFinite(u) ? regimeAt(u, umfE, ut) : "—";
}
function bindCalc() {
  $("calcParticle").addEventListener("change", () => {
    const p = PARTICLES[$("calcParticle").value];
    if (!p) return;
    $("calcRhos").value = p.rhos;
    $("calcDp").value = p.dpUm;
    updateCalc();
  });
  ["calcGas", "calcDp", "calcRhos", "calcEps", "calcU"].forEach((id) =>
    $(id).addEventListener("input", updateCalc));
  updateCalc();
}

// ---------- Tab 2: Ergun term explorer ----------
function drawErgunExplorer() {
  const cv = $("ergunCanvas");
  if (!cv || !cv.closest(".tab-panel").classList.contains("active")) return;
  const { ctx, w, h } = setupCanvas(cv);
  const u = parseFloat($("ergunU").value);
  const dp = parseFloat($("ergunDp").value) * 1e-6;
  const gas = GASES[$("ergunGas").value];
  const eps = parseFloat($("ergunEps").value);

  const uMax = 1.0;
  let yMax = 0;
  for (let i = 0; i <= 100; i++) {
    const uu = (i / 100) * uMax;
    yMax = Math.max(yMax, P.ergunDeltaPPerL(uu, gas.mu, gas.rho, eps, dp));
  }
  yMax *= 1.1;
  const L = 76, R = 18, T = 16, B = 46;
  const X = (uu) => L + (uu / uMax) * (w - L - R);
  const Y = (v) => h - B - (v / yMax) * (h - T - B);

  // Curves: viscous (blue), inertial (amber), total (red).
  const curves = [
    [(uu) => P.ergunViscousPerL(uu, gas.mu, eps, dp), "#4da3ff", "Viscous term"],
    [(uu) => P.ergunInertialPerL(uu, gas.rho, eps, dp), "#f5b301", "Inertial term"],
    [(uu) => P.ergunDeltaPPerL(uu, gas.mu, gas.rho, eps, dp), "#e00122", "Total"],
  ];
  curves.forEach(([fn, col]) => {
    ctx.strokeStyle = col;
    ctx.lineWidth = col === "#e00122" ? 2.5 : 2;
    ctx.beginPath();
    for (let i = 0; i <= 120; i++) {
      const uu = (i / 120) * uMax;
      if (i === 0) ctx.moveTo(X(uu), Y(fn(uu))); else ctx.lineTo(X(uu), Y(fn(uu)));
    }
    ctx.stroke();
  });

  // Marker at chosen u.
  ctx.strokeStyle = "#fff";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(X(u), Y(0));
  ctx.lineTo(X(u), Y(P.ergunDeltaPPerL(u, gas.mu, gas.rho, eps, dp)));
  ctx.stroke();

  // Legend.
  ctx.font = "12px sans-serif";
  ctx.textAlign = "left";
  curves.forEach(([fn, col, label], i) => {
    ctx.fillStyle = col;
    ctx.fillRect(L + 8, T + 8 + i * 18, 16, 4);
    ctx.fillStyle = "#e9e9e9";
    ctx.fillText(label, L + 30, T + 15 + i * 18);
  });

  const xTicks = [];
  for (let i = 0; i <= 5; i++) xTicks.push([X((i / 5) * uMax), ((i / 5) * uMax).toFixed(1)]);
  const yTicks = [];
  for (let i = 0; i <= 4; i++) {
    const v = (i / 4) * yMax;
    yTicks.push([Y(v), v >= 1000 ? (v / 1000).toFixed(1) + "k" : v.toFixed(0)]);
  }
  drawAxes(ctx, L, h - B, w - R, T, "Gas velocity  u  (m/s)", "ΔP/L  (Pa/m)", xTicks, yTicks);

  const visc = P.ergunViscousPerL(u, gas.mu, eps, dp);
  const iner = P.ergunInertialPerL(u, gas.rho, eps, dp);
  const tot = visc + iner;
  $("ergunVisc").textContent = fmt(visc, 0) + " Pa/m";
  $("ergunIner").textContent = fmt(iner, 0) + " Pa/m";
  $("ergunTot").textContent = fmt(tot, 0) + " Pa/m";
  $("ergunShare").textContent = fmt(100 * visc / tot, 0) + "% viscous";
}
function bindErgunExplorer() {
  const upd = () => {
    $("ergunUVal").textContent = fmt(parseFloat($("ergunU").value), 2) + " m/s";
    $("ergunDpVal").textContent = $("ergunDp").value + " μm";
    $("ergunGasVal").textContent = GASES[$("ergunGas").value].name;
    $("ergunEpsVal").textContent = fmt(parseFloat($("ergunEps").value), 2);
    drawErgunExplorer();
  };
  ["ergunU", "ergunDp", "ergunGas", "ergunEps"].forEach((id) =>
    $(id).addEventListener("input", upd));
  upd();
}

// ---------- Tab 4: Geldart chart ----------
const GEL_MEANINGS = [
  "C — Cohesive: interparticle forces dominate; the powder channels and plugs instead of fluidizing.",
  "A — Aeratable: fluidizes smoothly, expands well before bubbling. FCC catalyst lives here.",
  "B — Sand-like: bubbles as soon as it fluidizes; the classic bubbling bed.",
  "D — Spoutable: large and/or dense; needs high velocity and tends to spout.",
];
function drawGeldart() {
  const cv = $("geldartCanvas");
  if (!cv || !cv.closest(".tab-panel").classList.contains("active")) return;
  const { ctx, w, h } = setupCanvas(cv);
  const dpUm = parseFloat($("gelDp").value);
  const dRho = parseFloat($("gelDrho").value);

  // Log-log axes: dp 10..3000 um (x), dRho 100..10000 (y).
  const L = 76, R = 24, T = 20, B = 46;
  const X = (d) => L + (Math.log10(d) - 1) / (Math.log10(3000) - 1) * (w - L - R);
  const Y = (r) => h - B - (Math.log10(r) - 2) / (4 - 2) * (h - T - B);

  // Region shading.
  ctx.fillStyle = "rgba(77,163,255,0.06)";
  ctx.fillRect(L, T, w - L - R, h - B - T);
  // C region: dp < 30 um.
  ctx.fillStyle = "rgba(245,179,1,0.10)";
  ctx.fillRect(X(10), T, X(30) - X(10), h - B - T);
  // A/B and B/D boundary curves.
  const boundary = (k) => {
    ctx.strokeStyle = "#e9e9e9";
    ctx.lineWidth = 1.5;
    ctx.setLineDash([6, 4]);
    ctx.beginPath();
    for (let r = 100; r <= 10000; r *= 1.15) {
      const d = k * Math.sqrt(1000 / r);
      if (d < 10 || d > 3000) continue;
      const px = X(d), py = Y(r);
      if (r === 100) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();
    ctx.setLineDash([]);
  };
  boundary(100);   // A/B
  boundary(1000);  // B/D

  // Group labels.
  ctx.fillStyle = "#a8a8ad";
  ctx.font = "bold 15px sans-serif";
  ctx.textAlign = "center";
  ctx.fillText("C", X(18), Y(3000));
  ctx.fillText("A", X(55), Y(1200));
  ctx.fillText("B", X(400), Y(1200));
  ctx.fillText("D", X(1800), Y(1200));
  ctx.font = "11px sans-serif";
  ctx.fillText("A/B boundary (approx.)", X(120), Y(160));
  ctx.fillText("B/D boundary (approx.)", X(1200), Y(160));

  // User point.
  const group = P.geldartGroup(dpUm * 1e-6, dRho);
  ctx.fillStyle = "#e00122";
  ctx.beginPath();
  ctx.arc(X(dpUm), Y(dRho), 8, 0, 2 * Math.PI);
  ctx.fill();
  ctx.strokeStyle = "#fff";
  ctx.lineWidth = 2;
  ctx.stroke();

  const xTicks = [10, 30, 100, 300, 1000, 3000].map((d) => [X(d), String(d)]);
  const yTicks = [100, 1000, 10000].map((r) => [Y(r), String(r)]);
  drawAxes(ctx, L, h - B, w - R, T, "Particle diameter  d_p  (μm, log scale)", "ρ_s − ρ_g  (kg/m³, log scale)", xTicks, yTicks);

  $("gelGroup").textContent = "Group " + "CABD"[group];
  $("gelMeaning").textContent = GEL_MEANINGS[group];
}
function bindGeldart() {
  const upd = () => {
    $("gelDpVal").textContent = $("gelDp").value + " μm";
    $("gelDrhoVal").textContent = $("gelDrho").value + " kg/m³";
    drawGeldart();
  };
  ["gelDp", "gelDrho"].forEach((id) => $(id).addEventListener("input", upd));
  document.querySelectorAll("#gelPresets .chip").forEach((btn) => {
    btn.addEventListener("click", () => {
      const presets = {
        fcc: [75, 1500], sand: [500, 2600], fine: [20, 2600], coal: [1000, 1400],
      };
      const [d, r] = presets[btn.dataset.p];
      $("gelDp").value = d;
      $("gelDrho").value = r;
      document.querySelectorAll("#gelPresets .chip").forEach((c) => c.classList.toggle("active", c === btn));
      upd();
    });
  });
  upd();
}

// ---------- Tab 4: Terminal velocity explorer ----------
function drawUtExplorer() {
  const cv = $("utCanvas");
  if (!cv || !cv.closest(".tab-panel").classList.contains("active")) return;
  const { ctx, w, h } = setupCanvas(cv);
  const dp = parseFloat($("utDp").value) * 1e-6;
  const rhos = PARTICLES[$("utParticle").value].rhos;
  const gas = GASES[$("utGas").value];
  const ut = P.terminalVelocity(dp, rhos, gas.rho, gas.mu);
  const ReT = P.reynoldsTerminal(ut, dp, gas.rho, gas.mu);
  const umf = P.umfErgun(gas.mu, gas.rho, 0.4, dp, rhos, gas.rho);

  // Operating-window bar (log scale from 0.01 to 20 m/s).
  const L = 70, R = 20, y = h / 2 - 6;
  const X = (v) => L + (Math.log10(Math.max(v, 0.01)) + 2) / (Math.log10(20) + 2) * (w - L - R);
  ctx.fillStyle = "#222224";
  ctx.fillRect(L, y, w - L - R, 26);
  ctx.fillStyle = "rgba(55,201,120,0.35)";
  ctx.fillRect(X(umf), y, X(ut) - X(umf), 26);
  ctx.fillStyle = "#e00122";
  ctx.fillRect(X(umf) - 1.5, y - 6, 3, 38);
  ctx.fillStyle = "#4da3ff";
  ctx.fillRect(X(ut) - 1.5, y - 6, 3, 38);
  ctx.fillStyle = "#e9e9e9";
  ctx.font = "12px sans-serif";
  ctx.textAlign = "center";
  ctx.fillText("u_mf", X(umf), y - 12);
  ctx.fillText("u_t", X(ut), y - 12);
  ctx.fillStyle = "#37c978";
  ctx.fillText("operating window", (X(umf) + X(ut)) / 2, y + 44);
  ctx.fillStyle = "#a8a8ad";
  ctx.textAlign = "left";
  [0.01, 0.1, 1, 10].forEach((v) => ctx.fillText(v >= 1 ? v.toFixed(0) + " m/s" : v.toFixed(2) + " m/s", X(v) - 12, y + 62));

  const regime = ReT < 1 ? "Stokes" : ReT < 1000 ? "Intermediate" : "Newton";
  $("utVal").textContent = fmt(ut, 2) + " m/s";
  $("utRe").textContent = fmt(ReT, 1);
  $("utRegime").textContent = regime;
  $("utUmf").textContent = fmt(umf, 3) + " m/s";
  $("utWindow").textContent = "×" + fmt(ut / umf, 1);
}
function bindUtExplorer() {
  const upd = () => {
    $("utDpVal").textContent = $("utDp").value + " μm";
    drawUtExplorer();
  };
  ["utDp", "utParticle", "utGas"].forEach((id) => $(id).addEventListener("input", upd));
  upd();
}

// ---------- Tab 5: Bed expansion explorer ----------
function drawExpansion() {
  const cv = $("expandCanvas");
  if (!cv || !cv.closest(".tab-panel").classList.contains("active")) return;
  const { ctx, w, h } = setupCanvas(cv);
  const u = parseFloat($("expU").value);
  const dp = parseFloat($("expDp").value) * 1e-6;
  const rhos = PARTICLES[$("expParticle").value].rhos;
  const gas = GASES.air20;
  const epsMf = 0.4, Hmf = 0.5;
  const umf = P.umfErgun(gas.mu, gas.rho, epsMf, dp, rhos, gas.rho);
  const ut = P.terminalVelocity(dp, rhos, gas.rho, gas.mu);
  const ReT = P.reynoldsTerminal(ut, dp, gas.rho, gas.mu);

  let eps, H, state;
  if (u < umf) { eps = epsMf; H = Hmf; state = "Fixed bed"; }
  else if (u >= ut) { eps = 1; H = Hmf * 3; state = "Transport — bed emptying"; }
  else {
    eps = P.voidageAtVelocity(u, ut, ReT);
    H = P.bedHeight(Hmf, epsMf, eps);
    state = "Fluidized — expanded";
  }

  // Two columns: at u_mf and at u.
  const drawColumn = (x0, colW, HH, label) => {
    const base = h - 40;
    const scale = (h - 90) / (Hmf * 3);
    const top = base - HH * scale;
    ctx.strokeStyle = "#e9e9e9";
    ctx.lineWidth = 2;
    ctx.strokeRect(x0, top, colW, base - top);
    ctx.fillStyle = "rgba(224,1,34,0.55)";
    ctx.fillRect(x0 + 2, top + 2, colW - 4, base - top - 4);
    ctx.fillStyle = "#e9e9e9";
    ctx.font = "13px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(label, x0 + colW / 2, h - 18);
    ctx.fillText(fmt(HH, 2) + " m", x0 + colW / 2, top - 8);
  };
  drawColumn(w * 0.22, 110, Hmf, "At u_mf");
  drawColumn(w * 0.60, 110, H, "At u = " + fmt(u, 2) + " m/s");

  $("expUmf").textContent = fmt(umf, 3) + " m/s";
  $("expUt").textContent = fmt(ut, 2) + " m/s";
  $("expEps").textContent = fmt(eps, 3);
  $("expH").textContent = fmt(H, 2) + " m";
  $("expState").textContent = state;
}
function bindExpansion() {
  const upd = () => {
    $("expUVal").textContent = fmt(parseFloat($("expU").value), 2) + " m/s";
    $("expDpVal").textContent = $("expDp").value + " μm";
    drawExpansion();
  };
  ["expU", "expDp", "expParticle"].forEach((id) => $(id).addEventListener("input", upd));
  upd();
}

// ---------- Tab 6: FCC interactive diagram ----------
const FCC_INFO = {
  riser: "<b>Riser.</b> Heavy gas oil feed meets hot regenerated catalyst at the bottom. The mixture races upward together for only a few seconds — the vapour velocity is above the catalyst terminal velocity, so this is a transport (fast fluidization) regime. Cracking happens on the way up: big hydrocarbon molecules break into gasoline-range molecules on the zeolite catalyst.",
  reactor: "<b>Reactor / stripper.</b> At the top, cyclones separate the catalyst from the cracked vapours. Steam strips the last hydrocarbons off the spent catalyst before it leaves for the regenerator.",
  regen: "<b>Regenerator.</b> Air burns the coke off the spent catalyst at about 700 °C. The heat released here supplies the heat needed for cracking in the riser — the unit is thermally balanced. Flue gas leaves at the top; clean, hot catalyst returns to the riser.",
  frac: "<b>Fractionator.</b> The cracked vapours are condensed and separated into gasoline (the main product), LPG, light cycle oil (diesel-range), and heavy cycle oil.",
  catalyst: "<b>Catalyst circulation.</b> The catalyst loop is the heart of the unit: spent catalyst slides to the regenerator, clean catalyst flows back to the riser. The circulation rate controls the severity of cracking.",
  feed: "<b>Feed.</b> Heavy gas oil from crude distillation, preheated and atomised with steam, enters at the base of the riser.",
};
function bindFcc() {
  document.querySelectorAll("#fccSvg .fcc-part").forEach((el) => {
    el.addEventListener("click", () => {
      document.querySelectorAll("#fccSvg .fcc-part").forEach((e) => e.classList.remove("selected"));
      el.classList.add("selected");
      $("fccInfo").innerHTML = "<p>" + FCC_INFO[el.dataset.part] + "</p>";
    });
  });
}

// ---------- Draw everything (visible tab only) ----------
function drawAll() {
  drawDpGraph();
  drawBedFrame();
  drawErgunExplorer();
  drawGeldart();
  drawUtExplorer();
  drawExpansion();
}

// ---------- Init ----------
function initAll() {
  renderMath();
  // KaTeX loads with `defer`, so it arrives after this script runs.
  // Re-render once the document (and the deferred KaTeX scripts) are ready.
  document.addEventListener("DOMContentLoaded", renderMath);
  window.addEventListener("load", renderMath);
  initTabs();
  initParticles();
  bindDpControls();
  bindCalc();
  bindErgunExplorer();
  bindGeldart();
  bindUtExplorer();
  bindExpansion();
  bindFcc();
  $("animU").addEventListener("input", drawBedFrame);
  drawAll();
  bedLoop();
  window.addEventListener("resize", drawAll);
  if (!wasmOn) {
    console.warn("Running with the plain-JS physics fallback (WASM unavailable).");
  }
}

loadWasm();
