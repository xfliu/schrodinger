"""Assemble optimal_box_surface.csv and balance_law.json from the sweep output.

Reads the repaired dense-tier rows plus the LG-tier rows/detail, merges them into
one surface table, fits the two branches of the balance law, and picks the cells
worth an interval run.

The Moebius map is re-evaluated in post-processing wherever a different separator
rho is wanted: LG depends on (A0, A2, rho, c) in closed form and A0/A2 are
rho-independent, so no re-assembly is needed to change rho.
"""
import json, math
import numpy as np
import pandas as pd

CSHIFT = 1.0
DELTA_SEP = 0.0062257424      # published certified-mu2 separator deficit on Omega2 N=64
NUM = ["d_paper", "L", "Ly", "N", "Nprime", "D", "D_aux", "npanel", "t_star",
       "mu1N_Ritz", "eee2_Ritz", "lambda1D_Ritz", "mu2N_oee_Ritz", "eoe_Ritz",
       "rho_used", "LG_lower", "LG_lower_ideal_rho", "enclosure_width",
       "truncation_term", "LG_discretisation_term", "trunc_gap_discrete",
       "lambda_ref", "sigma_confinement", "max_a_over_L", "moment_tail_bound",
       "wall_seconds"]


def lg_mobius(A0, A2, rho, c=CSHIFT):
    """Single-test-vector Lehmann-Goerisch bound with A1 = 1 (float mirror of
    LGVerified.lg_bound): returns (lambda_lower, B, nu)."""
    rp = rho + c
    A = A0 - rp
    B = A0 - 2.0 * rp + rp * rp * A2
    nu = A / B
    return rp - rp / (1.0 - nu) - c, B, nu


def read_rows(path, broken_box=False):
    """Read a *_rows.csv.  `broken_box` repairs the pre-fix files whose unquoted
    box field split into three columns."""
    if not broken_box:
        return pd.read_csv(path)
    lines = open(path).read().strip().split("\n")
    hdr = lines[0].split(",")
    recs = []
    for ln in lines[1:]:
        f = ln.split(",")
        assert len(f) == len(hdr) + 2, (len(f), len(hdr))
        recs.append([f[0], f[1], ",".join(f[2:5])] + f[5:])
    return pd.DataFrame(recs, columns=hdr)


def coerce(df):
    for c in NUM:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["hypothesis_holds"] = df["hypothesis_holds"].astype(str).str.lower() == "true"
    df["key"] = df.system + " d=" + df.d_paper.map(lambda v: f"{v:g}")
    return df


def fit_power(x, y):
    """Least-squares fit of y = a * x^p (both positive); returns (p, a, r2)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return None
    lx, ly = np.log(x[m]), np.log(y[m])
    p, b = np.polyfit(lx, ly, 1)
    pred = p * lx + b
    ss = 1.0 - ((ly - pred) ** 2).sum() / max(((ly - ly.mean()) ** 2).sum(), 1e-300)
    return dict(exponent=float(p), prefactor=float(math.exp(b)), r2=float(ss), n=int(m.sum()))


def fit_exp(x, y):
    """Least-squares fit of y = a * exp(-k x); returns (k, a, r2)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = (y > 0) & np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return None
    k, b = np.polyfit(x[m], np.log(y[m]), 1)
    pred = k * x[m] + b
    ly = np.log(y[m])
    ss = 1.0 - ((ly - pred) ** 2).sum() / max(((ly - ly.mean()) ** 2).sum(), 1e-300)
    return dict(rate=float(-k), prefactor=float(math.exp(b)), r2=float(ss), n=int(m.sum()))
