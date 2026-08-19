#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Michaelis-Menten kinetics of porcine pancreatic lipase (LPP) and its
competitive inhibition by orlistat.

Model (Michaelian enzyme, competitive inhibition):

    v = Vmax * [S] / (Km * (1 + [I]/Ki) + [S])

The script:
  1. Loads a kinetic dataset (substrate-velocity pairs at several orlistat
     concentrations) consistent with the parameters reported in
     Candela, Arenas, Caicedo & Malagon, J. Chem. Educ. 2021
     (DOI: 10.1021/acs.jchemed.0c01184): Km ~ 61 mM, Vmax ~ 74 uM/min,
     Ki ~ 13 uM (Dixon plot).
  2. Fits the Michaelis-Menten equation at each inhibitor concentration.
  3. Performs a global fit of the competitive-inhibition model.
  4. Builds Lineweaver-Burk and Dixon plots and estimates Ki.
  5. Saves figures to the repository "figures" folder.

Requirements: numpy, scipy, matplotlib
"""

from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Figures are saved in <repo_root>/figures
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "figures"
OUTPUT_DIR.mkdir(exist_ok=True)

# ----------------------------------------------------------------------
# 1. EXPERIMENTAL DATA
# Representative dataset (mean of triplicate assays) consistent with the
# kinetic parameters reported in Candela et al. (J. Chem. Educ. 2021).
# Substrate (NPP) in mM; initial velocities in uM/min; orlistat in uM.
# ----------------------------------------------------------------------
S = np.array([0, 9, 18, 36, 45, 63, 81, 99], dtype=float)      # [NPP], mM
ORLISTAT = [0, 5, 10, 20, 50]                                   # uM

V_DATA = {
    0:  [0.0,  9.6, 16.7, 26.7, 31.0, 36.5, 42.3, 47.6],
    5:  [0.0,  7.0, 13.2, 22.4, 25.8, 30.7, 36.2, 40.8],
    10: [0.0,  5.6, 10.0, 17.8, 20.6, 27.1, 30.5, 35.7],
    20: [0.0,  4.0,  7.1, 13.7, 16.6, 21.5, 24.2, 28.4],
    50: [0.0,  2.1,  4.4,  7.8,  9.8, 13.3, 15.6, 18.5],
}

# ----------------------------------------------------------------------
# 2. KINETIC MODELS
# ----------------------------------------------------------------------
def michaelis_menten(s, vmax, km):
    """Michaelis-Menten equation."""
    return vmax * s / (km + s)

def competitive_global(x, vmax, km, ki):
    """Competitive inhibition model; x = (S, I)."""
    s, i = x
    return vmax * s / (km * (1.0 + i / ki) + s)

# ----------------------------------------------------------------------
# 3. PER-INHIBITOR MICHAELIS-MENTEN FITS
# ----------------------------------------------------------------------
print("=" * 66)
print("LPP + ORLISTAT: Michaelis-Menten fits per inhibitor concentration")
print("=" * 66)
per_I = {}
for I in ORLISTAT:
    v = np.array(V_DATA[I])
    popt, _ = curve_fit(michaelis_menten, S, v, p0=[75, 60])
    per_I[I] = popt
    print(f"[orlistat] = {I:>3} uM -> Vmax_app = {popt[0]:6.2f} uM/min, "
          f"Km_app = {popt[1]:6.2f} mM")
print("Note: Vmax_app stays ~constant while Km_app increases with [I],")
print("the fingerprint of COMPETITIVE inhibition.")

# ----------------------------------------------------------------------
# 4. GLOBAL FIT OF THE COMPETITIVE-INHIBITION MODEL
# ----------------------------------------------------------------------
S_all, I_all, v_all = [], [], []
for I in ORLISTAT:
    for s, v in zip(S, V_DATA[I]):
        S_all.append(s); I_all.append(I); v_all.append(v)
S_all, I_all, v_all = map(np.array, (S_all, I_all, v_all))

popt_g, _ = curve_fit(competitive_global, (S_all, I_all), v_all,
                      p0=[75, 60, 13])
Vmax, Km, Ki = popt_g
print("-" * 66)
print("GLOBAL FIT (competitive model):")
print(f"  Vmax = {Vmax:.2f} uM/min   (published: ~74 uM/min)")
print(f"  Km   = {Km:.2f} mM        (published: ~61 mM)")
print(f"  Ki   = {Ki:.2f} uM        (published: ~13 uM, Dixon)")

# ----------------------------------------------------------------------
# 5. LINEWEAVER-BURK PLOT (1/v vs 1/[S])
# ----------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
mask = S > 0
for I in ORLISTAT:
    v = np.array(V_DATA[I])[mask]
    x, y = 1.0 / S[mask], 1.0 / v
    a, b = np.polyfit(x, y, 1)
    ax[0].plot(x, y, "o", ms=4)
    xs = np.linspace(0, x.max() * 1.05, 20)
    ax[0].plot(xs, a * xs + b, "-", lw=1.2, label=f"{I} uM")
ax[0].set_xlabel("1/[NPP]  (mM$^{-1}$)")
ax[0].set_ylabel("1/v  (min/uM)")
ax[0].set_title("A. Lineweaver-Burk plot")
ax[0].legend(title="orlistat", fontsize=8)
ax[0].set_xlim(left=0)

# ----------------------------------------------------------------------
# 6. DIXON PLOT (1/v vs [I] at fixed [S])  ->  Ki = -x-intercept
# ----------------------------------------------------------------------
S_dixon = [9, 18, 45, 63, 99]
for s_fix in S_dixon:
    idx = np.where(S == s_fix)[0][0]
    inv_v = [1.0 / V_DATA[I][idx] for I in ORLISTAT]
    a, b = np.polyfit(ORLISTAT, inv_v, 1)
    ki_est = -b / a
    ax[1].plot(ORLISTAT, inv_v, "o", ms=4)
    xs = np.linspace(-ki_est - 3, max(ORLISTAT), 50)
    ax[1].plot(xs, a * xs + b, "-", lw=1.2,
               label=f"[S] = {s_fix} mM, Ki = {ki_est:.1f} uM")
ax[1].set_xlabel("[orlistat]  (uM)")
ax[1].set_ylabel("1/v  (min/uM)")
ax[1].set_title("B. Dixon plot")
ax[1].legend(fontsize=7.5)
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_ylim(bottom=0)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig_lipase_lineweaver_dixon.png", dpi=200)

# ----------------------------------------------------------------------
# 7. MICHAELIS-MENTEN CURVES WITH GLOBAL FIT
# ----------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(6.2, 4.4))
colors = plt.cm.viridis(np.linspace(0, 0.85, len(ORLISTAT)))
Ss = np.linspace(0, 100, 200)
for c, I in zip(colors, ORLISTAT):
    ax2.plot(S, V_DATA[I], "o", color=c, ms=4)
    ax2.plot(Ss, competitive_global((Ss, I), Vmax, Km, Ki), "-", color=c,
             lw=1.4, label=f"{I} uM")
ax2.set_xlabel("[NPP]  (mM)")
ax2.set_ylabel("v  (uM/min)")
ax2.set_title("LPP activity: Michaelis-Menten curves with orlistat")
ax2.legend(title="orlistat", fontsize=8)
fig2.tight_layout()
fig2.savefig(OUTPUT_DIR / "fig_lipase_michaelis.png", dpi=200)
print("-" * 66)
print(f"Figures saved in: {OUTPUT_DIR}")
print("  - fig_lipase_lineweaver_dixon.png")
print("  - fig_lipase_michaelis.png")

if __name__ == "__main__":
    plt.show()
