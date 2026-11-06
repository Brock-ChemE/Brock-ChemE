// Fluidized-bed physics core — compiled to WebAssembly with AssemblyScript.
// All functions use SI units. g = 9.81 m/s^2.
// Exports are plain f64 in / f64 out so the JS UI can call them directly.

// Ergun equation: pressure drop per unit bed length (Pa/m).
//   dP/L = 150*mu*(1-eps)^2*u / (eps^3*dp^2)   [viscous term]
//        + 1.75*(1-eps)*rho*u^2 / (eps^3*dp)    [inertial term]
export function ergunDeltaPPerL(u: f64, mu: f64, rho: f64, eps: f64, dp: f64): f64 {
  const e3: f64 = eps * eps * eps;
  const viscous: f64 = 150.0 * mu * (1.0 - eps) * (1.0 - eps) * u / (e3 * dp * dp);
  const inertial: f64 = 1.75 * (1.0 - eps) * rho * u * u / (e3 * dp);
  return viscous + inertial;
}

// Viscous (Kozeny-Carman) part of the Ergun pressure drop per unit length.
export function ergunViscousPerL(u: f64, mu: f64, eps: f64, dp: f64): f64 {
  const e3: f64 = eps * eps * eps;
  return 150.0 * mu * (1.0 - eps) * (1.0 - eps) * u / (e3 * dp * dp);
}

// Inertial (Burke-Plummer) part of the Ergun pressure drop per unit length.
export function ergunInertialPerL(u: f64, rho: f64, eps: f64, dp: f64): f64 {
  const e3: f64 = eps * eps * eps;
  return 1.75 * (1.0 - eps) * rho * u * u / (e3 * dp);
}

// Pressure drop across the bed at incipient fluidization (Pa):
// the buoyant weight of the bed per unit cross-sectional area.
export function bedWeightDeltaP(eps: f64, rhos: f64, rhog: f64, H: f64): f64 {
  return (1.0 - eps) * (rhos - rhog) * 9.81 * H;
}

// Minimum fluidization velocity from the Ergun equation.
// Set Ergun dP equal to the bed-weight dP and solve
//   B*u^2 + A*u - C = 0
//   A = 150*mu*(1-eps)^2/(eps^3*dp^2)
//   B = 1.75*(1-eps)*rho/(eps^3*dp)
//   C = (1-eps)*(rhos-rhog)*g
// Positive root: u = (-A + sqrt(A^2 + 4BC)) / (2B).
export function umfErgun(mu: f64, rho: f64, eps: f64, dp: f64, rhos: f64, rhog: f64): f64 {
  const e3: f64 = eps * eps * eps;
  const A: f64 = 150.0 * mu * (1.0 - eps) * (1.0 - eps) / (e3 * dp * dp);
  const B: f64 = 1.75 * (1.0 - eps) * rho / (e3 * dp);
  const C: f64 = (1.0 - eps) * (rhos - rhog) * 9.81;
  const disc: f64 = Math.sqrt(A * A + 4.0 * B * C);
  return (-A + disc) / (2.0 * B);
}

// Minimum fluidization velocity, small-particle closed form
// (viscous term dominates, Re_p < ~20):
//   u_mf = dp^2*(rhos-rhog)*g*eps^3 / (150*mu*(1-eps))
export function umfSmallParticle(mu: f64, eps: f64, dp: f64, rhos: f64, rhog: f64): f64 {
  return dp * dp * (rhos - rhog) * 9.81 * eps * eps * eps / (150.0 * mu * (1.0 - eps));
}

// Archimedes number: Ar = dp^3 * rho_g * (rho_s - rho_g) * g / mu^2
export function archimedes(dp: f64, rhos: f64, rhog: f64, mu: f64): f64 {
  return dp * dp * dp * rhog * (rhos - rhog) * 9.81 / (mu * mu);
}

// Wen-Yu correlation for minimum fluidization velocity (empirical cross-check):
//   Re_mf = sqrt(33.7^2 + 0.0408*Ar) - 33.7 ;  u_mf = Re_mf * mu / (dp * rho_g)
export function umfWenYu(dp: f64, rhos: f64, rhog: f64, mu: f64): f64 {
  const Ar: f64 = archimedes(dp, rhos, rhog, mu);
  const ReMf: f64 = Math.sqrt(33.7 * 33.7 + 0.0408 * Ar) - 33.7;
  return ReMf * mu / (dp * rhog);
}

// Drag coefficient from the Schiller-Naumann correlation.
export function dragCoefficient(Re: f64): f64 {
  if (Re < 1000.0) {
    return (24.0 / Re) * (1.0 + 0.15 * Math.pow(Re, 0.687));
  }
  return 0.44;
}

// Single-particle terminal velocity (m/s) by fixed-point iteration on
//   drag force = buoyant weight
//   u_t^2 = 4*dp*(rhos-rhog)*g / (3*rho_g*Cd),  Cd = Cd(Re), Re = rho_g*u_t*dp/mu
export function terminalVelocity(dp: f64, rhos: f64, rhog: f64, mu: f64): f64 {
  // Stokes estimate as the starting guess.
  let ut: f64 = dp * dp * (rhos - rhog) * 9.81 / (18.0 * mu);
  for (let i: i32 = 0; i < 200; i++) {
    const Re: f64 = rhog * ut * dp / mu;
    const Cd: f64 = dragCoefficient(Re);
    const utNew: f64 = Math.sqrt(4.0 * dp * (rhos - rhog) * 9.81 / (3.0 * rhog * Cd));
    if (Math.abs(utNew - ut) < 1e-9 * utNew) {
      ut = utNew;
      break;
    }
    ut = utNew;
  }
  return ut;
}

// Particle Reynolds number at the terminal velocity.
export function reynoldsTerminal(ut: f64, dp: f64, rhog: f64, mu: f64): f64 {
  return rhog * ut * dp / mu;
}

// Richardson-Zaki exponent n as a function of the terminal Reynolds number.
export function richardsonZakiN(ReT: f64): f64 {
  if (ReT < 0.2) return 4.65;
  if (ReT < 1.0) return 4.4 * Math.pow(ReT, -0.03);
  if (ReT < 500.0) return 4.4 * Math.pow(ReT, -0.1);
  return 2.4;
}

// Bed voidage at superficial velocity u from Richardson-Zaki: u/u_t = eps^n.
export function voidageAtVelocity(u: f64, ut: f64, ReT: f64): f64 {
  const n: f64 = richardsonZakiN(ReT);
  const ratio: f64 = u / ut;
  if (ratio <= 0.0) return 0.4;
  if (ratio >= 1.0) return 1.0;
  return Math.pow(ratio, 1.0 / n);
}

// Expanded bed height H = H_mf * (1 - eps_mf) / (1 - eps).
export function bedHeight(Hmf: f64, epsMf: f64, eps: f64): f64 {
  return Hmf * (1.0 - epsMf) / (1.0 - eps);
}

// Geldart group from particle diameter (m) and density difference (kg/m^3).
// Returns 0 = C (cohesive), 1 = A (aeratable), 2 = B (sand-like), 3 = D (spoutable).
// Boundaries are approximate (Geldart, 1973):
//   Group C: dp < 30 um
//   A/B boundary: dp_AB(um) = 100 * sqrt(1000 / dRho)
//   B/D boundary: dp_BD(um) = 1000 * sqrt(1000 / dRho)
export function geldartGroup(dp: f64, dRho: f64): i32 {
  const dpUm: f64 = dp * 1e6;
  if (dpUm < 30.0) return 0;
  const rhoRef: f64 = dRho > 1.0 ? dRho : 1.0;
  const dpAB: f64 = 100.0 * Math.sqrt(1000.0 / rhoRef);
  const dpBD: f64 = 1000.0 * Math.sqrt(1000.0 / rhoRef);
  if (dpUm < dpAB) return 1;
  if (dpUm < dpBD) return 2;
  return 3;
}
