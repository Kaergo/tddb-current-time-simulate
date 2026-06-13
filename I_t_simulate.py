#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fit SiC MOS gate-current I-t data with only FN tunneling and a defect
distribution.

The fitted current is:

    I(t) = K_FN * E_eff(t)^2 * exp(-B_FN / E_eff(t))

where the effective oxide field is reduced by occupied defect families:

    sigma_trap(t) = sum_k sigma_trap_k * (1 - exp(-c_k * t))
    E_eff(t) = E_0 - sigma_trap(t) / eps_ox

No TAT/PF component, TDDB generation term, or empirical comparison model is
used in the main fitting path.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re
from typing import Any

import matplotlib
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.integrate import quad

matplotlib.use("Agg")
import matplotlib.pyplot as plt


Q = 1.602176634e-19
HBAR = 1.054571817e-34
M0 = 9.1093837015e-31
PI = math.pi
EPS0 = 8.8541878128e-12

PHI_B_EV = 2.7
M_OX_REL = 0.42
NC_CM3 = 1.8e19
ND_CM3 = 1.0e16
M_SIC_REL = 0.29

def fermi_dirac_half(eta: float) -> float:
    eta = float(eta)
    if eta < -8.0:
        return float(np.exp(eta))
    if eta > 8.0:
        return float((4.0 / (3.0 * np.sqrt(np.pi))) * eta**1.5 * (1.0 + np.pi**2 / (8.0 * eta**2)))

    def integrand(eps):
        arg = np.clip(eps - eta, -100.0, 100.0)
        return np.sqrt(eps) / (1.0 + np.exp(arg))

    upper = max(80.0, eta + 80.0)
    value, _ = quad(integrand, 0.0, upper, limit=300)
    return float((2.0 / np.sqrt(np.pi)) * value)


def ec_bulk_minus_ef_eV(NC_cm3: float, ND_cm3: float, temperature_K: float) -> float:
    return float((8.617333262e-5 * temperature_K) * np.log(NC_cm3 / ND_cm3))


def surface_ns_cm3(
    NC_cm3: float,
    ND_cm3: float,
    temperature_K: float,
    phi_s_V: float,
) -> float:
    eta = (-ec_bulk_minus_ef_eV(NC_cm3, ND_cm3, temperature_K) + phi_s_V) / (8.617333262e-5 * temperature_K)
    return float(NC_cm3 * fermi_dirac_half(eta))


def get_exact_phi_s(Vg: float, Vfb: float, tox_nm: float, temperature_K: float) -> float:
    """Calculates surface potential using an internal Vg-phi_s lookup table without inverse solving."""
    kT_eV = 8.617333262e-5 * temperature_K
    ec_bulk = ec_bulk_minus_ef_eV(NC_CM3, ND_CM3, temperature_K)
    eps_sic = 10.32 * EPS0
    eps_ox = 3.9 * EPS0
    Cox = eps_ox / (tox_nm * 1e-9)
    
    def charge_integrand(phi):
        eta = (-ec_bulk + phi) / kT_eV
        n = NC_CM3 * fermi_dirac_half(eta)
        return n - ND_CM3
    
    phi_arr = np.linspace(0, 1.0, 50)
    Vg_arr = np.zeros(50)
    Vg_arr[0] = Vfb
    
    for i in range(1, len(phi_arr)):
        phi = phi_arr[i]
        integral, _ = quad(charge_integrand, 0, phi, limit=200)
        # qs in C/m^2 (integral is in V*cm^-3, so multiply by 1e6 for m^-3)
        qs = np.sqrt(2.0 * Q * eps_sic * integral * 1e6)
        Vg_arr[i] = Vfb + phi + qs / Cox
        
    return float(np.interp(Vg, Vg_arr, phi_arr))


def thermal_velocity_cm_s(temperature_K: float, m_eff_rel: float = M_SIC_REL) -> float:
    m = m_eff_rel * 9.1093837015e-31
    vth = math.sqrt(8.0 * 1.380649e-23 * temperature_K / (math.pi * m))
    return vth * 100.0  # m/s to cm/s



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit SiC MOS gate-current I-t data using only FN tunneling and defect distribution."
    )
    parser.add_argument("input_path", nargs="?", help="Optional input CSV or Excel file exported from the tester.")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Force demo mode even if a CSV path is provided.",
    )
    parser.add_argument("--ef-surf-eV", type=float, default=0.0, help="Fermi level relative to conduction band at surface")
    parser.add_argument("--tox-nm", type=float, default=50.0, help="Oxide thickness in nm.")
    parser.add_argument(
        "--sheet",
        type=str,
        default=None,
        help="Excel sheet name when the input file is .xlsx/.xls.",
    )
    parser.add_argument(
        "--time-column",
        type=str,
        default=None,
        help="Time column for Excel input, for example G or CH04_Time_s.",
    )
    parser.add_argument(
        "--current-column",
        type=str,
        default=None,
        help="Current column for Excel input, for example H or CH04_I_A.",
    )
    parser.add_argument(
        "--field-mvcm",
        type=float,
        default=8.3,
        help="Nominal oxide field in MV/cm. Used when --vg is omitted.",
    )
    parser.add_argument(
        "--vg",
        type=float,
        default=None,
        help="Stress voltage in V. If omitted, derived from field and tox.",
    )
    parser.add_argument("--width-um", type=float, default=200.0, help="Gate width in um.")
    parser.add_argument("--length-um", type=float, default=200.0, help="Gate length in um.")
    parser.add_argument(
        "--defect-families",
        type=int,
        default=4,
        help="Number of defect-distribution families used in sigma_trap(t).",
    )
    parser.add_argument(
        "--temperature-c",
        type=float,
        default=25.0,
        help="Device temperature in Celsius.",
    )
    parser.add_argument(
        "--vfb",
        type=float,
        default=0.0,
        help="Flatband voltage in V. Used to calculate oxide field if --vg is provided.",
    )
    parser.add_argument(
        "--scheme3",
        action="store_true",
        help="Enable Scheme 3 (discretized spatial exponential decay constraint) for fitting.",
    )
    parser.add_argument(
        "--scheme3-species",
        type=int,
        default=1,
        help="Number of spatial defect species for Scheme 3 (1, 2, 3, etc.)",
    )

    parser.add_argument(
        "--fit-bins",
        type=int,
        default=700,
        help="Number of log-time bins used for fitting; raw points are still used for reported metrics.",
    )
    parser.add_argument(
        "--eps-ox-r",
        type=float,
        default=3.9,
        help="Relative permittivity of SiO2 used in E_eff = E0 - sigma_trap/eps_ox.",
    )
    parser.add_argument(
        "--b-fn-vm",
        type=float,
        default=None,
        help="FN exponential coefficient B_FN in V/m. Default is derived from phi_b and m_ox.",
    )
    parser.add_argument(
        "--min-field-mvcm",
        type=float,
        default=0.1,
        help="Lower clipping value for E_eff during fitting, in MV/cm.",
    )
    parser.add_argument(
        "--fit-max-nfev",
        type=int,
        default=80000,
        help="Maximum function evaluations for each nonlinear least-squares start.",
    )
    parser.add_argument(
        "--ignore-first-point",
        action="store_true",
        help="Ignore the first data point before analysis.",
    )
    parser.add_argument(
        "--breakdown-ratio",
        type=float,
        default=100.0,
        help="Breakdown is detected when current exceeds this multiple of the recent median.",
    )
    parser.add_argument(
        "--breakdown-abs-a",
        type=float,
        default=1.0e-5,
        help="Absolute current threshold used together with --breakdown-ratio.",
    )
    parser.add_argument(
        "--breakdown-time-s",
        type=float,
        default=None,
        help="Manual breakdown onset time in seconds. If set, it overrides automatic detection.",
    )
    parser.add_argument(
        "--prefer-workbook-breakdown",
        action="store_true",
        help="For Excel workbooks, prefer the Breakdown sheet over automatic first-jump detection.",
    )
    parser.add_argument(
        "--demo-points",
        type=int,
        default=320,
        help="Number of time points in demo mode.",
    )
    parser.add_argument(
        "--demo-breakdown-s",
        type=float,
        default=1.0e4,
        help="Breakdown onset time in demo mode.",
    )
    parser.add_argument(
        "--demo-end-s",
        type=float,
        default=2.0e4,
        help="Last time point in demo mode.",
    )
    parser.add_argument(
        "--demo-noise-sigma",
        type=float,
        default=0.04,
        help="Log-normal relative noise used in demo mode.",
    )
    parser.add_argument("--output-dir", type=str, default=".", help="Directory for outputs.")
    parser.add_argument("--custom-bins-100", action="store_true", help="Use custom 10+90 binning scheme for fitting.")
    return parser.parse_args()


def derive_vg_and_field(tox_nm: float, field_mvcm: float, vg: float | None, vfb: float = 0.0) -> tuple[float, float]:
    if vg is None:
        vg_resolved = 0.1 * field_mvcm * tox_nm + vfb
        field_resolved = field_mvcm
    else:
        vg_resolved = vg
        field_resolved = 10.0 * (vg_resolved - vfb) / tox_nm
    return vg_resolved, field_resolved


def gate_area_cm2(width_um: float, length_um: float) -> float:
    return width_um * 1.0e-4 * length_um * 1.0e-4


def try_float(text: str) -> float | None:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def load_current_trace_csv(csv_path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    rows: list[list[str]] = []
    with csv_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        rows = list(csv.reader(f))

    times: list[float] = []
    currents: list[float] = []
    source_rows: list[int] = []
    for raw_idx, row in enumerate(rows, start=1):
        if row and row[0].strip() == "DataValue" and len(row) >= 3:
            t_val = try_float(row[1].strip())
            i_val = try_float(row[2].strip())
            if t_val is not None and i_val is not None:
                times.append(t_val)
                currents.append(i_val)
                source_rows.append(raw_idx)

    if not times:
        for raw_idx, row in enumerate(rows, start=1):
            numeric = [try_float(cell.strip()) for cell in row]
            numeric = [val for val in numeric if val is not None]
            if len(numeric) >= 2:
                times.append(numeric[0])
                currents.append(numeric[1])
                source_rows.append(raw_idx)

    if not times:
        raise ValueError(f"Could not parse time/current data from {csv_path}.")

    t_arr = np.asarray(times, dtype=float)
    i_arr = np.asarray(currents, dtype=float)
    row_arr = np.asarray(source_rows, dtype=int)
    mask = np.isfinite(t_arr) & np.isfinite(i_arr) & (t_arr >= 0.0) & (i_arr > 0.0)
    if np.count_nonzero(mask) < 3:
        raise ValueError(f"Not enough valid positive data points in {csv_path}.")

    t_arr = t_arr[mask]
    i_arr = i_arr[mask]
    row_arr = row_arr[mask]
    order = np.argsort(t_arr)
    return t_arr[order], i_arr[order], {
        "source_mode": "csv",
        "channel": None,
        "time_column_name": None,
        "current_column_name": None,
        "source_rows_1based": row_arr[order],
    }


def excel_column_index(spec: str) -> int | None:
    token = spec.strip().upper()
    if not token.isalpha():
        return None
    value = 0
    for ch in token:
        value = value * 26 + (ord(ch) - ord("A") + 1)
    return value - 1


def resolve_excel_column(columns: list[Any], spec: str | None, default_index: int) -> tuple[Any, int]:
    if spec is None:
        return columns[default_index], default_index

    idx = excel_column_index(spec)
    if idx is not None:
        if idx < 0 or idx >= len(columns):
            raise ValueError(f"Excel column {spec} is out of range.")
        return columns[idx], idx

    if spec in columns:
        return spec, columns.index(spec)

    for idx, name in enumerate(columns):
        if str(name).strip() == spec.strip():
            return name, idx
    raise ValueError(f"Cannot resolve Excel column spec {spec!r}.")


def infer_channel_from_name(name: Any) -> int | None:
    match = re.search(r"CH\s*0*(\d+)", str(name), flags=re.IGNORECASE)
    return int(match.group(1)) if match else None


def load_current_trace_excel(
    excel_path: Path,
    sheet: str | None,
    time_column: str | None,
    current_column: str | None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    df = pd.read_excel(excel_path, sheet_name=sheet if sheet is not None else 0)
    columns = list(df.columns)
    if len(columns) < 2:
        raise ValueError(f"Excel sheet {sheet!r} in {excel_path} has fewer than two columns.")

    time_col_name, time_idx = resolve_excel_column(columns, time_column, 0)
    current_col_name, current_idx = resolve_excel_column(columns, current_column, 1)
    if time_idx == current_idx:
        raise ValueError("Time column and current column must be different.")

    sub = df[[time_col_name, current_col_name]].copy()
    sub.columns = ["time", "current"]
    sub["source_row_1based"] = np.arange(len(df), dtype=int) + 2
    sub["time"] = pd.to_numeric(sub["time"], errors="coerce")
    sub["current"] = pd.to_numeric(sub["current"], errors="coerce")
    sub = sub.dropna()
    sub = sub[(sub["time"] >= 0.0) & (sub["current"] > 0.0)]
    if len(sub) < 3:
        raise ValueError(f"Not enough valid positive points in {excel_path} sheet {sheet!r}.")

    sub = sub.sort_values("time")
    channel = infer_channel_from_name(current_col_name) or infer_channel_from_name(time_col_name)
    return (
        sub["time"].to_numpy(dtype=float),
        sub["current"].to_numpy(dtype=float),
        {
            "source_mode": "excel",
            "sheet": sheet if sheet is not None else str(df.__class__.__name__),
            "channel": channel,
            "time_column_name": str(time_col_name),
            "current_column_name": str(current_col_name),
            "source_rows_1based": sub["source_row_1based"].to_numpy(dtype=int),
        },
    )


def load_current_trace(
    path: Path | str,
    sheet: str | None,
    time_column: str | None,
    current_column: str | None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return load_current_trace_csv(path)
    if suffix in {".xlsx", ".xls"}:
        return load_current_trace_excel(path, sheet=sheet, time_column=time_column, current_column=current_column)
    raise ValueError(f"Unsupported input file type: {path.suffix}")


def workbook_breakdown_time(excel_path: Path, channel: int | None) -> float | None:
    if channel is None:
        return None
    try:
        df = pd.read_excel(excel_path, sheet_name="Breakdown")
    except Exception:
        return None
    if "Channel" not in df.columns or "BreakdownTime_s" not in df.columns:
        return None
    row = df[pd.to_numeric(df["Channel"], errors="coerce") == float(channel)]
    if row.empty:
        return None
    value = pd.to_numeric(row.iloc[0]["BreakdownTime_s"], errors="coerce")
    return float(value) if pd.notna(value) else None


def generate_demo_trace(
    field_mvcm: float,
    area_cm2: float,
    defect_families: int,
    eps_ox_r: float,
    b_fn_v_m: float,
    n_points: int,
    breakdown_s: float,
    end_s: float,
    noise_sigma: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    t = np.geomspace(1.0e-3, end_s, n_points)
    ncomp = max(1, int(defect_families))
    e0_v_m = field_mvcm * 1.0e8
    eps_ox = eps_ox_r * EPS0
    min_field_v_m = 0.1e8
    tau_grid = np.geomspace(max(10.0, end_s / 2.0e3), max(100.0, end_s * 0.6), ncomp)
    rates = 1.0 / tau_grid
    weights = np.linspace(0.7, 1.3, ncomp)
    weights = weights / np.sum(weights)
    q_total = 0.13 * eps_ox * e0_v_m
    sigma_sat = q_total * weights
    sigma_trap = (1.0 - np.exp(-np.outer(t, rates))).dot(sigma_sat)
    e_eff = np.maximum(e0_v_m - sigma_trap / eps_ox, min_field_v_m)

    target_i0 = 8.0e-6 * max(area_cm2 / gate_area_cm2(200.0, 200.0), 1.0e-6)
    log_k = math.log(target_i0) - 2.0 * math.log(e0_v_m) + b_fn_v_m / e0_v_m
    current = np.exp(np.clip(log_k + 2.0 * np.log(e_eff) - b_fn_v_m / e_eff, -745, 709))

    rng = np.random.default_rng(20260406)
    noise = np.exp(rng.normal(loc=0.0, scale=noise_sigma, size=t.size))
    current = current * noise

    bd_idx = int(np.searchsorted(t, breakdown_s))
    if bd_idx < t.size:
        post_level = max(current[max(bd_idx - 1, 0)] * 2.0e4, 2.0e-4)
        post_noise = np.exp(rng.normal(loc=0.0, scale=0.20, size=t.size - bd_idx))
        current[bd_idx:] = post_level * post_noise

    return t, current, {
        "demo_K_FN_eff": float(math.exp(log_k)),
        "demo_E0_V_m": float(e0_v_m),
        "demo_B_FN_V_m": float(b_fn_v_m),
        "demo_Qdef_final_C_m2": float(sigma_trap[min(bd_idx, sigma_trap.size - 1)]),
        "demo_Eeff_final_V_m": float(e_eff[min(bd_idx, e_eff.size - 1)]),
    }


def recent_median(values: np.ndarray, end_idx: int, window: int) -> float:
    start = max(0, end_idx - window)
    return float(np.median(values[start:end_idx]))


def find_breakdown_index(
    current: np.ndarray,
    ratio: float,
    abs_current_a: float,
    window: int = 21,
    warmup: int = 50,
) -> int:
    if current.size <= warmup:
        return current.size
    for idx in range(warmup, current.size):
        ref = recent_median(current, idx, window)
        if ref <= 0.0 or not np.isfinite(ref):
            continue
        if current[idx] > max(abs_current_a, ratio * ref):
            return idx
    return current.size


def r2_log(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    a = np.log10(np.maximum(y_true, 1.0e-30))
    b = np.log10(np.maximum(y_pred, 1.0e-30))
    denom = np.sum((a - np.mean(a)) ** 2)
    if denom <= 0.0:
        return float("nan")
    return float(1.0 - np.sum((a - b) ** 2) / denom)


def cumulative_trapezoid_manual(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    if y.size != x.size:
        raise ValueError("y and x must have the same size.")
    out = np.zeros_like(y, dtype=float)
    if y.size < 2:
        return out
    dx = np.diff(x)
    avg = 0.5 * (y[:-1] + y[1:])
    out[1:] = np.cumsum(avg * dx)
    return out


def default_b_fn_vm(phi_b_ev: float = PHI_B_EV, m_ox_rel: float = M_OX_REL) -> float:
    phi_b_j = phi_b_ev * Q
    m_ox = m_ox_rel * M0
    coeff = (4.0 * math.sqrt(2.0 * m_ox)) / (3.0 * Q * HBAR)
    return coeff * (phi_b_j**1.5)


def barrier_width_nm(field_mvcm: float) -> float:
    field_vm = field_mvcm * 1.0e8
    phi_b_j = PHI_B_EV * Q
    return 1.0e9 * phi_b_j / (Q * field_vm)


def log_time_bin(t: np.ndarray, current: np.ndarray, bins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = np.isfinite(t) & np.isfinite(current) & (t > 0.0) & (current > 0.0)
    tb_src = t[mask]
    ib_src = current[mask]
    if tb_src.size < 3:
        raise ValueError("Need at least three positive-time positive-current points for log-time fitting.")

    order = np.argsort(tb_src)
    tb_src = tb_src[order]
    ib_src = ib_src[order]
    n_bins = max(3, min(int(bins), tb_src.size))
    t_min = float(tb_src[0])
    t_max = float(tb_src[-1])
    if not np.isfinite(t_min) or not np.isfinite(t_max) or t_max <= t_min:
        return tb_src, ib_src, np.ones_like(tb_src, dtype=int)

    pts_per_segment = 100
    t_1 = [t for t in tb_src if t <= 1000]
    edges_1 = np.geomspace(t_min, max(t_1) * (1.0 + 1e-12), pts_per_segment + 1) if t_1 else []
        
    t_2 = [t for t in tb_src if 1000 < t <= 10000]
    edges_2 = np.geomspace(max(1000, t_min), max(t_2) * (1.0 + 1e-12), pts_per_segment + 1) if t_2 else []
        
    t_3 = [t for t in tb_src if t > 10000]
    edges_3 = np.geomspace(max(10000, t_min), t_max * (1.0 + 1e-12), pts_per_segment + 1) if t_3 else []
        
    # Combine edges safely
    edges = []
    if len(edges_1) > 0:
        edges.extend(edges_1[:-1])
    if len(edges_2) > 0:
        edges.extend(edges_2[:-1])
    if len(edges_3) > 0:
        edges.extend(edges_3)
        
    edges = np.array(edges)
    n_bins_actual = len(edges) - 1

    which = np.searchsorted(edges, tb_src, side="right") - 1
    which = np.clip(which, 0, n_bins_actual - 1)

    t_out: list[float] = []
    i_out: list[float] = []
    counts: list[int] = []
    for idx in range(n_bins_actual):
        mask_idx = which == idx
        if np.any(mask_idx):
            t_out.append(float(np.median(tb_src[mask_idx])))
            i_out.append(float(np.median(ib_src[mask_idx])))
            counts.append(int(np.sum(mask_idx)))
    return np.asarray(t_out), np.asarray(i_out), np.asarray(counts, dtype=int)


def log_current_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(y_true) & np.isfinite(y_pred) & (y_true > 0.0) & (y_pred > 0.0)
    if np.count_nonzero(mask) < 3:
        return {
            "R2_log": float("nan"),
            "RMSE_ln": float("nan"),
            "MAE_ln": float("nan"),
            "typical_error_factor": float("nan"),
        }
    yt = np.log(y_true[mask])
    yp = np.log(y_pred[mask])
    residual = yt - yp
    denom = np.sum((yt - np.mean(yt)) ** 2)
    rmse = float(np.sqrt(np.mean(residual**2)))
    return {
        "R2_log": float(1.0 - np.sum(residual**2) / denom) if denom > 0.0 else float("nan"),
        "RMSE_ln": rmse,
        "MAE_ln": float(np.mean(np.abs(residual))),
        "typical_error_factor": float(np.exp(rmse)),
    }


def fn_defect_components(
    p: np.ndarray,
    t: np.ndarray,
    ncomp: int,
    e0_v_m: float,
    eps_ox: float,
    b_fn_v_m: float,
    min_field_v_m: float,
    scheme3: bool = False,
    ns: float = 1e20,
    vth: float = 1e7,
    temperature_K: float = 423.15,
    ef_surf_eV: float = 0.0,
    scheme3_species: int = 1,
    tox_nm: float = 45.0,
) -> dict[str, np.ndarray | float]:
    log_k = float(p[0])
    if scheme3:
        # WKB scheme parameters: [log_k0, log_Nt0, log_xd, log_tau00]
        log_Nt0 = float(p[1])
        log_xd = float(p[2])
        log_tau00 = float(p[3])

        Nt0_cm3 = np.exp(log_Nt0)
        Nt0_m3 = Nt0_cm3 * 1e6
        xd_m = np.exp(log_xd) * 1e-9
        tau00 = np.exp(log_tau00)
        tox_m = tox_nm * 1e-9

        # Grid setup
        n_grid = 80
        xmax_m = min(tox_m, 3.0e-9)
        dx_m = xmax_m / n_grid
        
        x_grid = np.linspace(dx_m/2, xmax_m - dx_m/2, n_grid)
        N_sheet = Nt0_m3 * np.exp(-x_grid / xd_m) * dx_m
        
        m_ox = M_OX_REL * M0
        wkb_coeff = 2.0 * np.sqrt(2.0 * m_ox * Q) / HBAR

        # ODE Solver setup
        num_steps = 1000
        t_grid = np.logspace(min(-4, np.log10(max(1e-6, t[0]))), np.log10(max(1e-3, t[-1])), num_steps)
        if t_grid[0] > 0:
            t_grid = np.insert(t_grid, 0, 0.0)

        f_trap = np.zeros_like(x_grid)
        sigma_trap_history = [0.0]
        S_FN_history = []
        
        # Initial WKB FN action calculation
        e_inj_0 = e0_v_m
        e_ox_local_v_m_0 = np.full_like(x_grid, e_inj_0)
        uc_local_0 = PHI_B_EV - np.cumsum(e_ox_local_v_m_0 * dx_m)
        S_FN_0 = 0.0
        for i, uc in enumerate(uc_local_0):
            if uc > 0:
                S_FN_0 += wkb_coeff * np.sqrt(uc) * dx_m
            else:
                break
        S_FN_history.append(S_FN_0)
        I0 = np.exp(log_k)

        for i in range(1, len(t_grid)):
            dt = t_grid[i] - t_grid[i-1]
            
            # Electrostatics
            sigma_i = Q * N_sheet * f_trap
            sigma_total = np.sum(sigma_i)
            shielding = np.sum((sigma_i / eps_ox) * ((tox_m - x_grid) / tox_m))
            e_inj = max(e0_v_m - shielding, min_field_v_m)
            
            e_ox_local_v_m = e_inj + np.cumsum(np.insert(sigma_i[:-1], 0, 0.0)) / eps_ox
            uc_local = PHI_B_EV - np.cumsum(e_ox_local_v_m * dx_m)
            
            # Update occupancies
            S_cap = wkb_coeff * np.cumsum(np.sqrt(np.maximum(0.0, uc_local)) * dx_m)
            tau_i = tau00 * np.exp(np.clip(S_cap, 0, 150))
            f_trap = 1.0 - (1.0 - f_trap) * np.exp(-dt / np.maximum(tau_i, 1e-30))
            
            # Record FN tunneling Action
            S_FN = 0.0
            for j, uc in enumerate(uc_local):
                if uc > 0:
                    S_FN += wkb_coeff * np.sqrt(uc) * dx_m
                else:
                    break
            
            sigma_trap_history.append(sigma_total)
            S_FN_history.append(S_FN)

        sigma_trap_history = np.array(sigma_trap_history)
        S_FN_history = np.array(S_FN_history)
        
        S_FN_interp = np.interp(t, t_grid, S_FN_history)
        current = I0 * np.exp(-np.clip(S_FN_interp - S_FN_0, -700, 700))
        sigma_trap = np.interp(t, t_grid, sigma_trap_history)
        e_eff = e0_v_m * np.ones_like(t) # Dummy fallback
        rates = np.zeros(1)
        sigma_sat = np.array([np.sum(Q * N_sheet)])
    else:
        sigma_sat = np.exp(p[1 : 1 + ncomp])
        rates = np.exp(p[1 + ncomp : 1 + 2 * ncomp])
        fill = 1.0 - np.exp(-np.outer(t, rates))
        sigma_trap = fill.dot(sigma_sat)
        e_eff = np.maximum(e0_v_m - sigma_trap / eps_ox, min_field_v_m)
        current = np.exp(np.clip(log_k + 2.0 * np.log(e_eff) - b_fn_v_m / e_eff, -745, 709))

    return {
        "current": current,
        "qdef_c_m2": sigma_trap,
        "e_eff_v_m": e_eff,
        "qcaps_c_m2": sigma_sat,
        "rates_s_inv": rates,
        "i_no_defect_a": current[0] if isinstance(current, np.ndarray) and len(current) > 0 else (I0 if scheme3 else np.exp(np.clip(log_k + 2.0 * np.log(e0_v_m) - b_fn_v_m / e0_v_m, -745, 709))),
        "k_fn_eff": np.exp(log_k)
    }


def fit_fn_defect_model(
    t_fit: np.ndarray,
    current_fit: np.ndarray,
    defect_families: int,
    e0_v_m: float,
    eps_ox: float,
    b_fn_v_m: float,
    min_field_v_m: float,
    max_nfev: int,
    scheme3: bool = False,
    ns: float = 1e20,
    vth: float = 1e7,
    temperature_K: float = 423.15,
    ef_surf_eV: float = 0.0,
    scheme3_species: int = 1,
    tox_nm: float = 45.0,
) -> dict[str, Any]:
    ncomp = max(1, int(defect_families))
    if e0_v_m <= min_field_v_m:
        raise ValueError("Nominal oxide field must be larger than --min-field-mvcm.")

    mask = np.isfinite(t_fit) & np.isfinite(current_fit) & (t_fit >= 0.0) & (current_fit > 0.0)
    t_local = t_fit[mask]
    i_local = current_fit[mask]
    if t_local.size < 3:
        raise ValueError("Not enough positive-current fit points.")

    order = np.argsort(t_local)
    t_local = t_local[order]
    i_local = i_local[order]
    log_i = np.log(i_local)

    first_i = float(i_local[0])
    log_k0 = math.log(first_i) - 2.0 * math.log(e0_v_m) + b_fn_v_m / e0_v_m
    q_upper = max(eps_ox * (e0_v_m - min_field_v_m) * 0.98, 1.0e-10)
    q_lower = 1.0e-14
    rate_lower = 1.0e-12
    rate_upper = 1.0e2

    positive_t = t_local[t_local > 0.0]
    t_min = float(np.min(positive_t)) if positive_t.size else 1.0e-6
    t_max = float(np.max(t_local))
    tau_min = max(t_min, 1.0e-6)
    tau_max = max(t_max, tau_min * 10.0)

    p0s: list[np.ndarray] = []
    if scheme3:
        # p = [log_k0, log_Nt0, log_xd, log_tau00]
        base_seed = [log_k0, np.log(1e19), np.log(0.5), np.log(1e-12)]
        p0s.append(np.array(base_seed))
        p0s.append(np.array([log_k0, np.log(1e18), np.log(1.0), np.log(1e-10)]))
        p0s.append(np.array([log_k0, np.log(1e20), np.log(0.2), np.log(1e-15)]))
        lb = np.array([log_k0 - 1.5, np.log(1e15), np.log(0.02), np.log(1e-30)])
        ub = np.array([log_k0 + 1.5, np.log(1e22), np.log(3.0), np.log(1e5)])
    else:
        base_tau = np.geomspace(tau_min, tau_max, ncomp)
        seed_rates = 1.0 / base_tau
        seed_q = np.full(ncomp, 1.0e-3)
        seed_arr = np.r_[log_k0, np.log(seed_q), np.log(seed_rates)]
        p0s.append(seed_arr)

        alt_rates = seed_rates * 2.0
        alt_q = np.full(ncomp, 5.0e-4)
        p0s.append(np.r_[log_k0, np.log(alt_q), np.log(alt_rates)])

        lb = np.r_[log_k0 - 30.0, np.log(np.full(ncomp, q_lower)), np.log(np.full(ncomp, rate_lower))]
        ub = np.r_[log_k0 + 30.0, np.log(np.full(ncomp, q_upper)), np.log(np.full(ncomp, rate_upper))]

    def residual(p: np.ndarray) -> np.ndarray:
        pred = fn_defect_components(
            p, t_local, ncomp, e0_v_m, eps_ox, b_fn_v_m, min_field_v_m,
            scheme3=scheme3, ns=ns, vth=vth, temperature_K=temperature_K, ef_surf_eV=ef_surf_eV,
            scheme3_species=scheme3_species, tox_nm=tox_nm
        )["current"]
        
        i_meas = np.exp(log_i)
        res = pred - i_meas
        
        i_range = np.max(i_meas) - np.min(i_meas)
        if i_range > 0:
            res = res / i_range
        
        # Apply physical weighting: early stages have tiny absolute drop but are physically crucial
        weight = np.ones_like(t_local)
        weight[t_local <= 1000] = 10.0
        weight[(t_local > 1000) & (t_local <= 10000)] = 3.0
        
        return res * weight

    best = None
    errors: list[str] = []
    for seed in p0s:
        try:
            res = least_squares(
                residual,
                seed,
                bounds=(lb, ub),
                max_nfev=max_nfev,
                x_scale="jac",
                verbose=2,
            )
        except Exception as exc:
            errors.append(str(exc))
            continue
        if best is None or res.cost < best.cost:
            best = res
            
        if scheme3:
            log_Nt0_res = res.x[1]
            log_xd_res = res.x[2]
            log_tau00_res = res.x[3]
            print(f"Seed: Nt0=10^{seed[1]/np.log(10):.1f}, x_decay=10^{seed[2]/np.log(10):.2f}, tau00=10^{seed[3]/np.log(10):.1f} -> "
                  f"Result: Nt0=10^{log_Nt0_res/np.log(10):.2f}, x_decay=10^{log_xd_res/np.log(10):.2f}, tau00=10^{log_tau00_res/np.log(10):.1f}, Cost={res.cost:.4e}")

    if best is None:
        joined = "; ".join(errors[-3:]) if errors else "unknown optimizer failure"
        raise RuntimeError(f"FN-defect fit failed: {joined}")

    return {
        "params": best.x,
        "cost": float(best.cost),
        "status": int(best.status),
        "message": str(best.message),
        "nfev": int(best.nfev),
        "ncomp": ncomp,
        "e0_v_m": float(e0_v_m),
        "eps_ox_f_m": float(eps_ox),
        "b_fn_v_m": float(b_fn_v_m),
        "min_field_v_m": float(min_field_v_m),
    }


def family_rows(sigma_sat: np.ndarray, rates: np.ndarray) -> list[dict[str, float | str]]:
    order = np.argsort(rates)[::-1]
    rows: list[dict[str, float | str]] = []
    for rank, idx in enumerate(order, start=1):
        qcap = float(sigma_sat[idx])
        rate = float(rates[idx])
        rows.append(
            {
                "family": f"G{rank}",
                "Qcap_C_m2": qcap,
                "c_s_inv": rate,
                "tau_s": float(1.0 / rate) if rate > 0 else float("inf"),
                "Ncap_cm_minus2": float(qcap / Q / 1.0e4),
            }
        )
    return rows


def json_sanitize(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): json_sanitize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_sanitize(v) for v in value]
    return value


def analyze_trace(
    t_pre: np.ndarray,
    i_pre: np.ndarray,
    width_um: float,
    length_um: float,
    tox_nm: float,
    field_mvcm: float,
    defect_families: int,
    fit_bins: int,
    eps_ox_r: float,
    b_fn_v_m: float,
    min_field_mvcm: float,
    max_nfev: int,
    scheme3: bool = False,
    temperature_c: float = 25.0,
    vg: float = 0.0,
    vfb: float = 0.0,
    scheme3_species: int = 1,
) -> dict[str, Any]:
    area_cm2 = gate_area_cm2(width_um, length_um)
    temperature_K = temperature_c + 273.15
    
    # EXACT surface potential via lookup table / integration
    phi_s_V = get_exact_phi_s(vg, vfb, tox_nm, temperature_K)
    ef_surf_eV = phi_s_V - ec_bulk_minus_ef_eV(NC_CM3, ND_CM3, temperature_K)
    ns = surface_ns_cm3(NC_CM3, ND_CM3, temperature_K, phi_s_V)
    
    # Correct oxide field calculation to account for flatband and surface potential
    if vg > 0.0:
        e0_v_m = 1.0e9 * (vg - vfb - phi_s_V) / tox_nm
    else:
        e0_v_m = field_mvcm * 1.0e8
        
    eps_ox = eps_ox_r * EPS0
    min_field_v_m = min_field_mvcm * 1.0e8
    tb, ib, bin_count = log_time_bin(t_pre, i_pre, fit_bins)
    
    vth = thermal_velocity_cm_s(temperature_K)

    fit = fit_fn_defect_model(
        t_fit=tb,
        current_fit=ib,
        defect_families=defect_families,
        e0_v_m=e0_v_m,
        eps_ox=eps_ox,
        b_fn_v_m=b_fn_v_m,
        min_field_v_m=min_field_v_m,
        max_nfev=max_nfev,
        scheme3=scheme3,
        ns=ns,
        vth=vth,
        temperature_K=temperature_K,
        ef_surf_eV=ef_surf_eV,
        scheme3_species=scheme3_species,
        tox_nm=tox_nm,
    )

    ncomp = fit["ncomp"]
    raw = fn_defect_components(fit["params"], t_pre, ncomp, e0_v_m, eps_ox, b_fn_v_m, min_field_v_m, scheme3=scheme3, ns=ns, vth=vth, temperature_K=temperature_K, ef_surf_eV=ef_surf_eV, scheme3_species=scheme3_species, tox_nm=tox_nm)
    binned = fn_defect_components(fit["params"], tb, ncomp, e0_v_m, eps_ox, b_fn_v_m, min_field_v_m, scheme3=scheme3, ns=ns, vth=vth, temperature_K=temperature_K, ef_surf_eV=ef_surf_eV, scheme3_species=scheme3_species, tox_nm=tox_nm)
    sigma_sat = np.asarray(raw["qcaps_c_m2"], dtype=float)
    rates = np.asarray(raw["rates_s_inv"], dtype=float)
    families = family_rows(sigma_sat, rates)

    model_current = np.asarray(raw["current"], dtype=float)
    no_defect_current = np.full_like(t_pre, float(raw["i_no_defect_a"]), dtype=float)
    qinj_meas_c = cumulative_trapezoid_manual(i_pre, t_pre)
    qinj_model_c = cumulative_trapezoid_manual(model_current, t_pre)
    qinj_no_defect_c = cumulative_trapezoid_manual(no_defect_current, t_pre)

    raw_metrics = log_current_metrics(i_pre, model_current)
    binned_metrics = log_current_metrics(ib, np.asarray(binned["current"], dtype=float))
    qdef_end = float(np.asarray(raw["qdef_c_m2"])[-1])
    e_end = float(np.asarray(raw["e_eff_v_m"])[-1])

    return {
        "area_cm2": area_cm2,
        "field_mvcm": field_mvcm,
        "barrier_width_nm": barrier_width_nm(field_mvcm),
        "tox_nm": tox_nm,
        "eps_ox_r": eps_ox_r,
        "eps_ox_f_m": eps_ox,
        "b_fn_v_m": b_fn_v_m,
        "e0_v_m": e0_v_m,
        "min_field_v_m": min_field_v_m,
        "t_pre": t_pre,
        "i_pre": i_pre,
        "j_pre": i_pre / area_cm2,
        "t_bin": tb,
        "i_bin": ib,
        "bin_count": bin_count,
        "i_model": model_current,
        "j_model": model_current / area_cm2,
        "i_no_defect": no_defect_current,
        "j_no_defect": no_defect_current / area_cm2,
        "qdef_c_m2": np.asarray(raw["qdef_c_m2"], dtype=float),
        "e_eff_v_m": np.asarray(raw["e_eff_v_m"], dtype=float),
        "e_eff_mv_cm": np.asarray(raw["e_eff_v_m"], dtype=float) / 1.0e8,
        "qinj_meas_c": qinj_meas_c,
        "qinj_model_c": qinj_model_c,
        "qinj_no_defect_c": qinj_no_defect_c,
        "qinj_meas_c_cm2": qinj_meas_c / area_cm2,
        "qinj_model_c_cm2": qinj_model_c / area_cm2,
        "qinj_no_defect_c_cm2": qinj_no_defect_c / area_cm2,
        "qbd_meas_c": float(qinj_meas_c[-1]),
        "qbd_model_c": float(qinj_model_c[-1]),
        "qbd_no_defect_c": float(qinj_no_defect_c[-1]),
        "qbd_meas_c_cm2": float(qinj_meas_c[-1] / area_cm2),
        "qbd_model_c_cm2": float(qinj_model_c[-1] / area_cm2),
        "qbd_no_defect_c_cm2": float(qinj_no_defect_c[-1] / area_cm2),
        "fit": fit,
        "k_fn_eff": float(raw["k_fn_eff"]),
        "i_no_defect_a": float(raw["i_no_defect_a"]),
        "families": families,
        "qcap_total_c_m2": float(np.sum(sigma_sat)),
        "ncap_total_cm_minus2": float(np.sum(sigma_sat) / Q / 1.0e4),
        "qdef_end_c_m2": qdef_end,
        "n_def_occupied_end_cm_minus2": float(qdef_end / Q / 1.0e4),
        "e_eff_end_v_m": e_end,
        "e_eff_end_mv_cm": e_end / 1.0e8,
        "field_drop_end_mv_cm": (e0_v_m - e_end) / 1.0e8,
        "raw_metrics": raw_metrics,
        "binned_metrics": binned_metrics,
        "r2_log": raw_metrics["R2_log"],
    }


def save_outputs(
    output_dir: Path,
    stem: str,
    analysis: dict[str, Any],
    breakdown_index: int,
    breakdown_time_s: float,
    breakdown_source: str,
    source_mode: str,
    vg_resolved: float,
    tox_nm: float,
    width_um: float,
    length_um: float,
    input_meta: dict[str, Any],
    scheme3: bool = False,
    scheme3_species: int = 1,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_m = analysis["raw_metrics"]
    bin_m = analysis["binned_metrics"]
    summary = f"""Input stem: {stem}
Source mode: {source_mode}

Model family:
- FN tunneling current only: I = K_FN * E_eff^2 * exp(-B_FN / E_eff)
- Defect families only enter by trapped charge screening the oxide field.
- sigma_trap(t) = sum_k sigma_trap_k * (1 - exp(-c_k * t))
- E_eff(t) = E0 - sigma_trap(t) / eps_ox
- No TAT, PF, TDDB generation, or empirical comparison model is included.

Geometry / stress:
tox = {tox_nm:.3f} nm
width = {width_um:.3f} um
length = {length_um:.3f} um
area = {analysis['area_cm2']:.6e} cm^2
stress voltage = {vg_resolved:.6f} V
oxide field = {analysis['field_mvcm']:.6f} MV/cm
total FN barrier width W = {analysis['barrier_width_nm']:.6f} nm
eps_ox_r = {analysis['eps_ox_r']:.6g}
B_FN = {analysis['b_fn_v_m']:.6e} V/m

Breakdown exclusion:
source = {breakdown_source}
index (0-based filtered trace) = {breakdown_index}
point (1-based filtered trace) = {breakdown_index + 1}
time = {breakdown_time_s:.6f} s
input columns = {input_meta.get('time_column_name')} / {input_meta.get('current_column_name')}
channel = {input_meta.get('channel')}
source row (1-based original file) = {input_meta.get('breakdown_source_row_1based')}

FN-defect fit:
optimizer status = {analysis['fit']['status']}
optimizer message = {analysis['fit']['message']}
function evaluations = {analysis['fit']['nfev']}
defect families = {analysis['fit']['ncomp']}
K_FN_eff = {analysis['k_fn_eff']:.6e}
FN no-defect current = {analysis['i_no_defect_a']:.6e} A
total Qcap = {analysis['qcap_total_c_m2']:.6e} C/m^2
total Ncap = {analysis['ncap_total_cm_minus2']:.6e} cm^-2
occupied sigmaTrap at final pre-breakdown point = {analysis['qdef_end_c_m2']:.6e} C/m^2
occupied Ndef at final pre-breakdown point = {analysis['n_def_occupied_end_cm_minus2']:.6e} cm^-2
E_eff final = {analysis['e_eff_end_mv_cm']:.6f} MV/cm
field drop final = {analysis['field_drop_end_mv_cm']:.6f} MV/cm

Fit quality in log(current):
binned R2_log = {bin_m['R2_log']:.6f}
binned RMSE_ln = {bin_m['RMSE_ln']:.6e}
binned typical error factor = {bin_m['typical_error_factor']:.6f}
raw R2_log = {raw_m['R2_log']:.6f}
raw RMSE_ln = {raw_m['RMSE_ln']:.6e}
raw typical error factor = {raw_m['typical_error_factor']:.6f}

Charge-to-breakdown integrated to the last pre-breakdown point:
Measured QBD total      = {analysis['qbd_meas_c']:.6e} C
Measured QBD density    = {analysis['qbd_meas_c_cm2']:.6e} C/cm^2
FN-defect model QBD         = {analysis['qbd_model_c']:.6e} C
FN-defect model QBD density = {analysis['qbd_model_c_cm2']:.6e} C/cm^2
FN no-defect QBD            = {analysis['qbd_no_defect_c']:.6e} C
FN no-defect QBD density    = {analysis['qbd_no_defect_c_cm2']:.6e} C/cm^2
"""
    family_lines = ["", "Defect family parameters:"]
    for i, row in enumerate(analysis["families"]):
        if scheme3 and i >= 5 and i < len(analysis["families"]) - 5:
            if i == 5:
                family_lines.append("... (intermediate spatial bins omitted) ...")
            continue
        family_lines.append(
            f"{row['family']}: Qcap={row['Qcap_C_m2']:.6e} C/m^2, "
            f"Ncap={row['Ncap_cm_minus2']:.6e} cm^-2, "
            f"c={row['c_s_inv']:.6e} s^-1, tau={row['tau_s']:.6e} s"
        )
    summary += "\n".join(family_lines) + "\n"
    if scheme3 and "params" in analysis["fit"]:
        p = analysis["fit"]["params"]
        log_Nt0_val = p[1]
        log_xd_val = p[2]
        log_tau00_val = p[3]
        summary += f"\nWKB Trap-Charging-Induced FN Current Relaxation Model parameters:\n"
        summary += f"  Peak Defect Density Nt0 = 10^{log_Nt0_val/np.log(10):.6f} cm^-3\n"
        summary += f"  Spatial Decay Length xd = 10^{log_xd_val/np.log(10):.6f} nm\n"
        summary += f"  Pre-exponential Time Constant tau00 = 10^{log_tau00_val/np.log(10):.6f} s\n"
    (output_dir / f"{stem}_summary.txt").write_text(summary, encoding="utf-8")

    data = {
        "time_s": analysis["t_pre"],
        "measured_current_A": analysis["i_pre"],
        "measured_current_density_A_per_cm2": analysis["j_pre"],
        "model_FN_defect_current_A": analysis["i_model"],
        "model_FN_no_defect_current_A": analysis["i_no_defect"],
        "model_FN_defect_J_A_per_cm2": analysis["j_model"],
        "model_FN_no_defect_J_A_per_cm2": analysis["j_no_defect"],
        "Qdef_C_m2": analysis["qdef_c_m2"],
        "Eeff_V_m": analysis["e_eff_v_m"],
        "Eeff_MV_cm": analysis["e_eff_mv_cm"],
        "cum_charge_measured_C": analysis["qinj_meas_c"],
        "cum_charge_measured_C_per_cm2": analysis["qinj_meas_c_cm2"],
        "cum_charge_FN_defect_C": analysis["qinj_model_c"],
        "cum_charge_FN_defect_C_per_cm2": analysis["qinj_model_c_cm2"],
        "cum_charge_FN_no_defect_C": analysis["qinj_no_defect_c"],
        "cum_charge_FN_no_defect_C_per_cm2": analysis["qinj_no_defect_c_cm2"],
    }
    np.savetxt(
        output_dir / f"{stem}_fit_data.csv",
        np.column_stack(list(data.values())),
        delimiter=",",
        header=",".join(data.keys()),
        comments="",
    )
    pd.DataFrame(analysis["families"]).to_csv(output_dir / f"{stem}_defect_families.csv", index=False)
    pd.DataFrame(
        {
            "time_s": analysis["t_bin"],
            "median_current_A": analysis["i_bin"],
            "bin_count": analysis["bin_count"],
        }
    ).to_csv(output_dir / f"{stem}_fit_bins.csv", index=False)
    fit_summary = {
        "input_stem": stem,
        "source_mode": source_mode,
        "model": "FN_tunneling_with_defect_charge_distribution",
        "geometry": {
            "tox_nm": tox_nm,
            "width_um": width_um,
            "length_um": length_um,
            "area_cm2": analysis["area_cm2"],
        },
        "stress": {
            "vg_V": vg_resolved,
            "field_MV_cm": analysis["field_mvcm"],
            "E0_V_m": analysis["e0_v_m"],
            "B_FN_V_m": analysis["b_fn_v_m"],
        },
        "breakdown": {
            "source": breakdown_source,
            "index_0based": breakdown_index,
            "time_s": breakdown_time_s,
            "source_row_1based": input_meta.get("breakdown_source_row_1based"),
        },
        "fit": {
            "optimizer_status": analysis["fit"]["status"],
            "optimizer_message": analysis["fit"]["message"],
            "nfev": analysis["fit"]["nfev"],
            "K_FN_eff": analysis["k_fn_eff"],
            "I_FN_no_defect_A": analysis["i_no_defect_a"],
            "raw_metrics": analysis["raw_metrics"],
            "binned_metrics": analysis["binned_metrics"],
            "families": analysis["families"],
        },
    }
    if scheme3 and "params" in analysis["fit"]:
        p = analysis["fit"]["params"]
        n_sp = scheme3_species
        fit_summary["fit"]["dynamic_coulomb_global_parameters"] = {
            "tau0_s": float(10**p[1]),
            "tau1_s": float(10**p[2]),
            "x_decay_nm": float(p[3]),
        }
        fit_summary["fit"]["wkb_fn_relaxation_parameters"] = {
            "log_k0": float(p[0]),
            "Nt0_cm3": float(math.exp(p[1])),
            "xd_nm": float(math.exp(p[2])),
            "tau00_s": float(math.exp(p[3]))
        }
    (output_dir / f"{stem}_fit_summary.json").write_text(
        json.dumps(json_sanitize(fit_summary), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    plt.figure(figsize=(8.0, 5.0))
    plt.semilogy(analysis["t_bin"], analysis["i_bin"], "o", markersize=3, label="Sampled Measured", color='blue', alpha=0.6)
    plt.semilogy(analysis["t_pre"], analysis["i_model"], linewidth=2, label="FN + defect distribution", color='red')
    plt.semilogy(analysis["t_pre"], analysis["i_no_defect"], "--", linewidth=1.5, label="FN no defect", color='green')
    plt.xlabel("Time (s)")
    plt.ylabel("Current (A)")
    plt.title("Pre-breakdown current: measured vs FN-defect model")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / f"{stem}_combined.png", dpi=220)
    plt.close()

    # Log-Log plot to inspect early stage fitting
    plt.figure(figsize=(8.0, 5.0))
    plt.loglog(analysis["t_bin"], analysis["i_bin"], "o", markersize=3, label="Sampled Measured", color='blue', alpha=0.6)
    plt.loglog(analysis["t_pre"], analysis["i_model"], linewidth=2, label="FN + defect distribution", color='red')
    plt.loglog(analysis["t_pre"], analysis["i_no_defect"], "--", linewidth=1.5, label="FN no defect", color='green')
    plt.xlabel("Time (s) - Log Scale")
    plt.ylabel("Current (A) - Log Scale")
    plt.title("Pre-breakdown current: early stage fitting (log-log)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / f"{stem}_combined_loglog.png", dpi=220)
    plt.close()

    def plot_interval(t_min, t_max, filename, title):
        if t_min == 0:
            mask_bin = (analysis["t_bin"] >= t_min) & (analysis["t_bin"] <= t_max)
        else:
            mask_bin = (analysis["t_bin"] > t_min) & (analysis["t_bin"] <= t_max)

        if not np.any(mask_bin):
            return 0.0, 0.0
            
        t_bin_mask = analysis["t_bin"][mask_bin]
        i_bin_mask = analysis["i_bin"][mask_bin]
        
        # Calculate model current at exact binned points
        log_i_mod = np.interp(t_bin_mask, analysis["t_pre"], np.log(np.maximum(analysis["i_model"], 1.0e-300)))
        i_model_bin = np.exp(log_i_mod)
        
        # Calculate NRMSE for the binned points in this interval in linear scale
        r2, rmse, nrmse = 0.0, 0.0, 0.0
        if len(t_bin_mask) > 2:
            i_meas = i_bin_mask
            i_mod = i_model_bin
            
            rmse = np.sqrt(np.mean((i_meas - i_mod)**2))
            
            i_range = np.max(i_meas) - np.min(i_meas)
            nrmse = rmse / i_range if i_range > 0 else 0.0
            
            text_str = f"$NRMSE$ = {nrmse:.2%}"
        else:
            text_str = ""
            
        plt.figure(figsize=(8.0, 5.0))
        plt.loglog(t_bin_mask, i_bin_mask, "o", markersize=4, label="Sampled Measured", color='blue', alpha=0.7)
        plt.loglog(t_bin_mask, i_model_bin, "r-", linewidth=2.5, label="FN + defect distribution")
        plt.xlabel("Time (s)")
        plt.ylabel("Current (A)")
        plt.title(title)
        
        if text_str:
            props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
            plt.gca().text(0.05, 0.05, text_str, transform=plt.gca().transAxes, fontsize=11,
                           verticalalignment='bottom', bbox=props)
                           
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig(output_dir / filename, dpi=220)
        plt.close()
        return nrmse

    nrmse_1 = plot_interval(0, 1000, f"{stem}_interval_1.png", "Fitting Interval 1: 0 to 1000 s")
    nrmse_2 = plot_interval(1000, 10000, f"{stem}_interval_2.png", "Fitting Interval 2: 1000 to 10000 s")
    nrmse_3 = plot_interval(10000, np.max(analysis["t_pre"]), f"{stem}_interval_3.png", "Fitting Interval 3: > 10000 s")

    avg_nrmse = (nrmse_1 + nrmse_2 + nrmse_3) / 3.0

    interval_info = f"\nInterval-specific fitting quality (NRMSE of linear current):\n"
    interval_info += f"  Interval 1 (0 - 1000 s): NRMSE = {nrmse_1:.2%}\n"
    interval_info += f"  Interval 2 (1000 - 10000 s): NRMSE = {nrmse_2:.2%}\n"
    interval_info += f"  Interval 3 (> 10000 s): NRMSE = {nrmse_3:.2%}\n"
    interval_info += f"  Average of three intervals: NRMSE = {avg_nrmse:.2%}\n"

    print(interval_info.strip())

    summary_path = output_dir / f"{stem}_summary.txt"
    try:
        existing_summary = summary_path.read_text(encoding="utf-8")
    except Exception:
        existing_summary = ""
    summary_path.write_text(existing_summary + "\n" + interval_info, encoding="utf-8")

    plt.figure(figsize=(8.0, 5.0))
    residual = np.log(np.maximum(analysis["i_pre"], 1.0e-300)) - np.log(
        np.maximum(analysis["i_model"], 1.0e-300)
    )
    plt.semilogx(analysis["t_pre"], residual, ".", markersize=2)
    plt.axhline(0.0, color="black", linewidth=1.0)
    plt.xlabel("Time (s)")
    plt.ylabel("ln(I_measured) - ln(I_model)")
    plt.title("FN-defect residual")
    plt.tight_layout()
    plt.savefig(output_dir / f"{stem}_residual.png", dpi=220)
    plt.close()

    plt.figure(figsize=(8.0, 5.0))
    plt.semilogx(analysis["t_pre"], analysis["qdef_c_m2"], linewidth=2)
    plt.xlabel("Time (s)")
    plt.ylabel("sigma_trap (C/m$^2$)")
    plt.title("Occupied defect charge")
    plt.tight_layout()
    plt.savefig(output_dir / f"{stem}_qdef.png", dpi=220)
    plt.close()

    plt.figure(figsize=(8.0, 5.0))
    plt.semilogx(analysis["t_pre"], analysis["e_eff_mv_cm"], linewidth=2)
    plt.xlabel("Time (s)")
    plt.ylabel("E_eff (MV/cm)")
    plt.title("Effective oxide field screened by defects")
    plt.tight_layout()
    plt.savefig(output_dir / f"{stem}_eeff.png", dpi=220)
    plt.close()

    plt.figure(figsize=(8.0, 5.0))
    plt.plot(analysis["t_pre"], analysis["qinj_meas_c_cm2"], label="Measured Q/A")
    plt.plot(analysis["t_pre"], analysis["qinj_model_c_cm2"], label="FN-defect model Q/A")
    plt.plot(analysis["t_pre"], analysis["qinj_no_defect_c_cm2"], label="FN no-defect Q/A")
    plt.xlabel("Time (s)")
    plt.ylabel("Injected charge density (C/cm$^2$)")
    plt.title("Cumulative injected charge before breakdown")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / f"{stem}_qbd.png", dpi=220)
    plt.close()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    area_cm2 = gate_area_cm2(args.width_um, args.length_um)
    b_fn_v_m = float(args.b_fn_vm) if args.b_fn_vm is not None else default_b_fn_vm()

    demo_mode = args.demo or args.input_path is None
    if demo_mode:
        vg_resolved, field_mvcm = derive_vg_and_field(args.tox_nm, args.field_mvcm, args.vg, args.vfb)
        t, i, demo_meta = generate_demo_trace(
            field_mvcm=field_mvcm,
            area_cm2=area_cm2,
            defect_families=args.defect_families,
            eps_ox_r=args.eps_ox_r,
            b_fn_v_m=b_fn_v_m,
            n_points=args.demo_points,
            breakdown_s=args.demo_breakdown_s,
            end_s=args.demo_end_s,
            noise_sigma=args.demo_noise_sigma,
        )
        stem = "demo_fn_defect_current"
        source_mode = "demo"
        input_meta = {
            "time_column_name": "demo_time_s",
            "current_column_name": "demo_current_A",
            "channel": None,
        }
        breakdown_time_hint = None
    else:
        input_path = Path(args.input_path)
        t, i, input_meta = load_current_trace(
            input_path,
            sheet=args.sheet,
            time_column=args.time_column,
            current_column=args.current_column,
        )
        vg_resolved, field_mvcm = derive_vg_and_field(args.tox_nm, args.field_mvcm, args.vg, args.vfb)
        stem = input_path.stem.replace(" ", "_").replace("/", "_")
        source_mode = input_meta["source_mode"]
        demo_meta = {}
        breakdown_time_hint = (
            workbook_breakdown_time(input_path, input_meta.get("channel"))
            if input_path.suffix.lower() in {".xlsx", ".xls"}
            else None
        )

    if args.ignore_first_point and t.size > 1:
        t = t[1:]
        i = i[1:]
        source_rows = input_meta.get("source_rows_1based")
        if isinstance(source_rows, np.ndarray) and source_rows.size == t.size + 1:
            input_meta["source_rows_1based"] = source_rows[1:]

    if args.breakdown_time_s is not None:
        breakdown_index = int(np.searchsorted(t, args.breakdown_time_s, side="left"))
        breakdown_time_s = float(args.breakdown_time_s)
        breakdown_source = "manual_breakdown_time"
    elif args.prefer_workbook_breakdown and breakdown_time_hint is not None:
        breakdown_index = int(np.searchsorted(t, breakdown_time_hint, side="left"))
        breakdown_time_s = float(breakdown_time_hint)
        breakdown_source = "workbook_breakdown_sheet"
    else:
        breakdown_index = find_breakdown_index(
            current=i,
            ratio=args.breakdown_ratio,
            abs_current_a=args.breakdown_abs_a,
        )
        if breakdown_index >= t.size and breakdown_time_hint is not None:
            breakdown_index = int(np.searchsorted(t, breakdown_time_hint, side="left"))
            breakdown_time_s = float(breakdown_time_hint)
            breakdown_source = "workbook_breakdown_sheet_fallback"
        else:
            breakdown_time_s = float(t[breakdown_index]) if breakdown_index < t.size else float("nan")
            breakdown_source = "auto_first_sharp_rise"
    t_pre = t[:breakdown_index]
    i_pre = i[:breakdown_index]
    source_rows = input_meta.get("source_rows_1based")
    if isinstance(source_rows, np.ndarray) and breakdown_index < source_rows.size:
        input_meta["breakdown_source_row_1based"] = int(source_rows[breakdown_index])
    else:
        input_meta["breakdown_source_row_1based"] = None
    if t_pre.size < 3:
        raise ValueError("Not enough pre-breakdown points remain for analysis.")

    analysis = analyze_trace(
        t_pre=t_pre,
        i_pre=i_pre,
        width_um=args.width_um,
        length_um=args.length_um,
        tox_nm=args.tox_nm,
        field_mvcm=field_mvcm,
        defect_families=args.defect_families,
        fit_bins=args.fit_bins,
        eps_ox_r=args.eps_ox_r,
        b_fn_v_m=b_fn_v_m,
        min_field_mvcm=args.min_field_mvcm,
        max_nfev=args.fit_max_nfev,
        scheme3=args.scheme3,
        temperature_c=args.temperature_c,
        vg=vg_resolved,
        vfb=args.vfb,
        scheme3_species=getattr(args, "scheme3_species", 1),
    )

    save_outputs(
        output_dir=output_dir,
        stem=stem,
        analysis=analysis,
        breakdown_index=breakdown_index,
        breakdown_time_s=breakdown_time_s,
        breakdown_source=breakdown_source,
        source_mode=source_mode,
        vg_resolved=vg_resolved,
        tox_nm=args.tox_nm,
        width_um=args.width_um,
        length_um=args.length_um,
        input_meta=input_meta,
        scheme3=args.scheme3,
        scheme3_species=getattr(args, "scheme3_species", 1),
    )

    print("Done.")
    print(f"Mode: {source_mode}")
    print(f"Output directory: {output_dir.resolve()}")
    print(f"oxide field = {field_mvcm:.4f} MV/cm")
    print(f"input columns = {input_meta.get('time_column_name')} / {input_meta.get('current_column_name')}")
    print(f"channel = {input_meta.get('channel')}")
    print(f"breakdown source = {breakdown_source}")
    print(f"breakdown index (0-based) = {breakdown_index}")
    print(f"breakdown point (1-based filtered trace) = {breakdown_index + 1}")
    print(f"breakdown time = {breakdown_time_s:.6f} s")
    print(f"breakdown source row (1-based original file) = {input_meta.get('breakdown_source_row_1based')}")
    print(f"barrier width W = {analysis['barrier_width_nm']:.4f} nm")
    print(f"B_FN = {analysis['b_fn_v_m']:.6e} V/m")
    print(f"K_FN_eff = {analysis['k_fn_eff']:.6e}")
    print(f"FN no-defect current = {analysis['i_no_defect_a']:.6e} A")
    print(f"defect families = {analysis['fit']['ncomp']}")
    print(f"total Ncap = {analysis['ncap_total_cm_minus2']:.6e} cm^-2")
    print(f"final E_eff = {analysis['e_eff_end_mv_cm']:.6f} MV/cm")
    print(f"field drop final = {analysis['field_drop_end_mv_cm']:.6f} MV/cm")
    print(f"raw R2_log = {analysis['raw_metrics']['R2_log']:.6f}")
    print(f"binned R2_log = {analysis['binned_metrics']['R2_log']:.6f}")
    print(f"Measured QBD total      = {analysis['qbd_meas_c']:.6e} C")
    print(f"Measured QBD density    = {analysis['qbd_meas_c_cm2']:.6e} C/cm^2")
    print(f"FN-defect model QBD     = {analysis['qbd_model_c']:.6e} C")
    print(f"FN-defect QBD density   = {analysis['qbd_model_c_cm2']:.6e} C/cm^2")
    if demo_meta:
        print(
            "Demo source: "
            f"K_FN={demo_meta['demo_K_FN_eff']:.6e}, "
            f"Qdef_final={demo_meta['demo_Qdef_final_C_m2']:.6e} C/m^2, "
            f"Eeff_final={demo_meta['demo_Eeff_final_V_m'] / 1.0e8:.6f} MV/cm"
        )


if __name__ == "__main__":
    main()
