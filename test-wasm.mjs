// Node verification: instantiate build/release.wasm and compare exported
// physics functions against the spec check values (sand, dp=500 um,
// rhos=2600, air 20 C rho=1.2, mu=1.81e-5, eps=0.4, H=0.5 m).
import fs from "node:fs";

const bytes = fs.readFileSync(new URL("./build/release.wasm", import.meta.url));
const { instance } = await WebAssembly.instantiate(bytes, {});
const p = instance.exports;

const dp = 500e-6, rhos = 2600, rhog = 1.2, mu = 1.81e-5, eps = 0.4, H = 0.5;

function check(name, actual, expected, tolPct) {
  const err = Math.abs(actual - expected) / expected * 100;
  const ok = err <= tolPct;
  console.log(`${ok ? "PASS" : "FAIL"}  ${name}: actual=${actual.toFixed(4)} expected~${expected} err=${err.toFixed(2)}%`);
  return ok;
}

let allOk = true;
allOk &= check("Ergun dP/L at u=0.1 m/s (Pa/m)", p.ergunDeltaPPerL(0.1, mu, rhog, eps, dp), 6469, 5);
allOk &= check("Bed-weight dP_mf (Pa)", p.bedWeightDeltaP(eps, rhos, rhog, H), 7648, 5);
allOk &= check("u_mf Ergun quadratic (m/s)", p.umfErgun(mu, rhog, eps, dp, rhos, rhog), 0.22, 5);
allOk &= check("u_mf Wen-Yu (m/s)", p.umfWenYu(dp, rhos, rhog, mu), 0.20, 5);
allOk &= check("u_mf small-particle closed form (m/s)", p.umfSmallParticle(mu, eps, dp, rhos, rhog), 0.250, 5);
allOk &= check("Terminal velocity (m/s)", p.terminalVelocity(dp, rhos, rhog, mu), 3.9, 5);
allOk &= check("Archimedes number", p.archimedes(dp, rhos, rhog, mu), 11672, 5);
allOk &= check("Richardson-Zaki n at Re_t", p.richardsonZakiN(p.reynoldsTerminal(p.terminalVelocity(dp, rhos, rhog, mu), dp, rhog, mu)), 2.71, 5);

// Geldart group checks (0=C, 1=A, 2=B, 3=D)
function checkGroup(name, actual, expected) {
  const ok = actual === expected;
  console.log(`${ok ? "PASS" : "FAIL"}  ${name}: actual=${actual} expected=${expected}`);
  return ok;
}
allOk &= checkGroup("Geldart FCC catalyst (75 um, dRho 1500) -> A", p.geldartGroup(75e-6, 1500), 1);
allOk &= checkGroup("Geldart sand (500 um, dRho 2600) -> B", p.geldartGroup(500e-6, 2600), 2);
allOk &= checkGroup("Geldart fine powder (20 um) -> C", p.geldartGroup(20e-6, 2600), 0);
allOk &= checkGroup("Geldart coal (1000 um, dRho 1400) -> D", p.geldartGroup(1000e-6, 1400), 3);

// Sanity: Ergun terms sum to total; voidage/height behave sensibly.
const visc = p.ergunViscousPerL(0.1, mu, eps, dp);
const iner = p.ergunInertialPerL(0.1, rhog, eps, dp);
console.log(`INFO  Ergun terms at u=0.1: viscous=${visc.toFixed(1)} inertial=${iner.toFixed(1)} sum=${(visc + iner).toFixed(1)}`);
const ut = p.terminalVelocity(dp, rhos, rhog, mu);
console.log(`INFO  Re_t=${p.reynoldsTerminal(ut, dp, rhog, mu).toFixed(1)} voidage at u=1 m/s=${p.voidageAtVelocity(1.0, ut, p.reynoldsTerminal(ut, dp, rhog, mu)).toFixed(3)} bed height at u=1 (H_mf=0.5)=${p.bedHeight(0.5, 0.4, p.voidageAtVelocity(1.0, ut, p.reynoldsTerminal(ut, dp, rhog, mu))).toFixed(3)} m`);

console.log(allOk ? "\nALL CHECKS PASSED" : "\nSOME CHECKS FAILED");
process.exit(allOk ? 0 : 1);
