// Energy-aligned oxide trap capture model core for SiC MOS I-t simulation
// This file is intentionally DOM-free so it can be imported into UI or tested independently.

export const Q = 1.602176634e-19;
export const EPS0 = 8.8541878128e-12;
export const HBAR = 1.054571817e-34;
export const M0 = 9.1093837015e-31;

export function clamp(x, lo, hi) {
  return Math.max(lo, Math.min(hi, x));
}

export function smoothStep(Eev, widthEv) {
  const w = Math.max(widthEv, 1e-6);
  const z = clamp(Eev / w, -60, 60);
  return 1 / (1 + Math.exp(-z));
}

export function buildTrapGrid(state) {
  const n = Math.max(8, Math.round(state.currentTrapGridN || 80));
  const toxM = state.toxNm * 1e-9;
  const xmaxM = Math.min(toxM, Math.max(0.05e-9, state.currentXmaxNm * 1e-9));
  const dxM = xmaxM / n;
  const Nt0M3 = Math.max(0, state.currentNt0Cm3) * 1e6;
  const xdM = Math.max(0.02e-9, state.currentXdecayNm * 1e-9);
  const traps = [];
  for (let i = 0; i < n; i += 1) {
    const xM = (i + 0.5) * dxM;
    const NtM3 = Nt0M3 * Math.exp(-xM / xdM);
    const NsheetM2 = NtM3 * dxM;
    traps.push({xM, xNm: xM * 1e9, dxM, NtM3, NsheetM2, f: 0, sigma: 0, EtRelEv: 0, Hcb: 0, tau: Infinity});
  }
  return traps;
}

export function solveOxideElectrostatics(state, traps) {
  const epsOx = state.epsRel * EPS0;
  const toxM = state.toxNm * 1e-9;
  let EinjVm = state.e0Vm;
  let sigmaTotal = 0;
  for (const tr of traps) {
    tr.sigma = Q * tr.NsheetM2 * tr.f;
    sigmaTotal += tr.sigma;
    EinjVm -= (tr.sigma / epsOx) * ((toxM - tr.xM) / Math.max(toxM, 1e-30));
  }
  return {EinjVm, eEffVm: EinjVm, sigmaTotal, traps};
}

export function EofX(state, electro, xM) {
  const epsOx = state.epsRel * EPS0;
  let E = electro.EinjVm;
  for (const tr of electro.traps) {
    if (tr.xM < xM) E += tr.sigma / epsOx;
  }
  return E;
}

export function barrierUcEv(state, electro, xM) {
  const epsOx = state.epsRel * EPS0;
  let dropV = electro.EinjVm * xM;
  for (const tr of electro.traps) {
    if (tr.xM < xM) dropV += (tr.sigma / epsOx) * (xM - tr.xM);
  }
  return state.phiB - dropV;
}

export function wkbActionToXAtEnergy(state, electro, xEndM, energyRefEv) {
  if (!(xEndM > 0)) return 0;
  const mEff = state.mOx * M0;
  const steps = Math.max(24, Math.ceil(xEndM / (0.02e-9)));
  const dx = xEndM / steps;
  let integral = 0;
  for (let i = 0; i < steps; i += 1) {
    const x = (i + 0.5) * dx;
    const U = Math.max(barrierUcEv(state, electro, x) - energyRefEv, 0);
    if (U > 0) integral += Math.sqrt(U * Q) * dx;
  }
  return (2 * Math.sqrt(2 * mEff) / HBAR) * integral;
}

export function wkbActionFn(state, electro) {
  return wkbActionToXAtEnergy(state, electro, state.toxNm * 1e-9, 0);
}

export function captureRateForTrap(state, electro, tr) {
  const UtrapEv = barrierUcEv(state, electro, tr.xM);
  const EtRelEv = UtrapEv - state.currentDeltaEtEv;
  const Hcb = smoothStep(EtRelEv, state.currentAlignWidthEv);

  const Sgap = wkbActionToXAtEnergy(state, electro, tr.xM, 0);
  const cGap = Math.exp(-clamp(Sgap, 0, 120)) / state.currentTau0Gap;

  let cCb = 0;
  let Scb = 0;
  if (Hcb > 1e-8 && EtRelEv > -10 * state.currentAlignWidthEv) {
    const eRef = Math.max(EtRelEv, 0);
    Scb = wkbActionToXAtEnergy(state, electro, tr.xM, eRef);
    cCb = Math.exp(-clamp(Scb, 0, 120)) / state.currentTau0Cb;
  }

  const cTotal = Hcb * cCb + (1 - Hcb) * cGap;
  tr.EtRelEv = EtRelEv;
  tr.Hcb = Hcb;
  tr.Sgap = Sgap;
  tr.Scb = Scb;
  tr.cCb = cCb;
  tr.cGap = cGap;
  tr.tau = cTotal > 0 ? 1 / cTotal : Infinity;
  return cTotal;
}

export function fnInitialCurrent(state) {
  const E = Math.max(state.e0Vm, 1);
  const mEff = state.mOx * M0;
  const phiJ = state.phiB * Q;
  const prefactor = (Q * Q * E * E) / (8 * Math.PI * HBAR * phiJ);
  const action = (4 * Math.sqrt(2 * mEff) * Math.pow(phiJ, 1.5)) / (3 * Q * HBAR * E);
  return state.areaM2 * prefactor * Math.exp(clamp(-action, -745, 80));
}

export function simulateCurrentRelaxation(state, times) {
  const traps = buildTrapGrid(state);
  const electro0 = solveOxideElectrostatics(state, traps);
  const S0 = wkbActionFn(state, electro0);
  const I0 = Number.isFinite(state.currentI0) && state.currentI0 > 0 ? state.currentI0 : fnInitialCurrent(state);
  const out = [];
  let tPrev = 0;

  for (const t of times) {
    let dtLeft = Math.max(0, t - tPrev);
    while (dtLeft > 0) {
      const electro = solveOxideElectrostatics(state, traps);
      let rateMax = 0;
      const rates = traps.map((tr) => {
        const c = captureRateForTrap(state, electro, tr);
        rateMax = Math.max(rateMax, c);
        return c;
      });
      const tauMin = rateMax > 0 ? 1 / rateMax : 1e99;
      const dt = Math.min(dtLeft, Math.max(dtLeft / 40, 0.25 * Math.max(tauMin, 1e-15)));
      if (!(dt > 0) || !Number.isFinite(dt)) break;
      for (let i = 0; i < traps.length; i += 1) {
        traps[i].f = 1 - (1 - traps[i].f) * Math.exp(-rates[i] * dt);
        traps[i].f = clamp(traps[i].f, 0, 1);
      }
      dtLeft -= dt;
    }
    const electro = solveOxideElectrostatics(state, traps);
    for (const tr of traps) captureRateForTrap(state, electro, tr);
    const S = wkbActionFn(state, electro);
    const I = I0 * Math.exp(clamp(-(S - S0), -745, 80));
    const totalNsheet = traps.reduce((s, tr) => s + tr.NsheetM2, 0);
    const avgHcb = traps.reduce((s, tr) => s + tr.Hcb * tr.NsheetM2, 0) / Math.max(totalNsheet, 1e-300);
    out.push({t, current: I, noDefect: I0, ratio: I / Math.max(I0, 1e-300), Sfn: S, eEffVm: electro.EinjVm, eEffMvCm: electro.EinjVm / 1e8, sigmaTrap: electro.sigmaTotal, avgHcb, traps: traps.map((tr) => ({...tr}))});
    tPrev = t;
  }
  return out;
}
