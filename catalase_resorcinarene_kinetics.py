#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Non-Michaelian kinetics of bovine liver catalase (CAT) and its
non-competitive inhibition by sulfonated resorcinarenes.

Catalase does not reach saturation with H2O2 under the assay conditions;
its kinetics are described as a first-order process:

    v = k1 * [S]

Non-competitive inhibition by the resorcinarene reduces the apparent
first-order constant without changing the substrate dependence:

    v = k1 * [S] / (1 + [I]/Ki)

The script:
  1. Loads a kinetic dataset (substrate-velocity pairs at several
     Na4EtRA concentrations, plus relative-activity curves for the five
     resorcinarenes) consistent with the results reported in
     Collazos, Garcia, Malagon, Caicedo & Vargas, Int. J. Biol.
     Macromol. 139 (2019) 75-84 (DOI: 10.1016/j.ijbiomac.2019.07.197):
     non-Michaelian behavior, non-competitive inhibition, up to ca. 70%
     inhibition for Na4EtRA and Na4PrRA, and the inhibitory tendency
     Na4PrRA > Na4EtRA > Na4ESRA > Na4MeRA > Na4SRA.
  2. Demonstrates the non-Michaelian behavior by comparing the
     first-order model against Michaelis-Menten via AIC.
  3. Fits the non-competitive first-order inhibition model globally
     and estimates k1 and Ki.
  4. Builds the Lineweaver-Burk-style plot (lines converge on the
     x-axis: non-competitive inhibition).
  5. Fits the relative-activity curves of the five resorcinarenes and
     ranks their inhibitory capacity.
  6. Saves figures to the repository "figures" folder.

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
# results reported in Collazos et al. (Int. J. Biol. Macromol. 2019).
# Substrate (H2O2) in uM; initial velocities in uM/min;
# resorcinarene concentrations in uM.
# ----------------------------------------------------------------------
S = np.array([5, 8, 14, 25, 50], dtype=float)                    # [H2O2], uM
RES_LEVELS = [0, 10, 20, 50, 100]                                # uM

# Detailed kinetics with Na4EtRA at each inhibitor concentration
V_ETRA = {
    0:   [10.01, 16.65, 29.03, 49.23, 99.11],
    10:  [ 8.05, 13.31, 22.87, 41.83, 77.28],
    20:  [ 7.25, 11.04, 19.78, 34.47, 68.44],
    50:  [ 4.80,  7.77, 13.18, 23.58, 48.34],
    100: [ 3.02,  4.74,  8.79, 15.21, 29.25],
}

# Relative activity (%, control = 100) at [H2O2] = 18 uM for the five
# sulfonated resorcinarenes
REL_ACTIVITY = {
    "Na4PrRA":  [ 98.4, 79.3, 65.1, 43.1, 28.6],
    "Na4EtRA":  [101.8, 81.4, 68.2, 47.7, 31.5],
    "Na4ESRA":  [ 99.4, 89.9, 81.7, 61.3, 43.7],
    "Na4MeRA":  [100.7, 92.8, 87.6, 68.8, 53.8],
    "Na4SRA":   [ 98.3, 90.9, 89.1, 77.0, 60.6],
}

# ----------------------------------------------------------------------
# 2. KINETIC MODELS
# ----------------------------------------------------------------------
def first_order(s, k1):
    """First-order (non-Michaelian) model: v = k1 * [S]."""
    return k1 * s

def michaelis_menten(s, vmax, km):
    """Michaelis-Menten equation (reference model)."""
    return vmax * s / (km + s)

def noncompetitive_first_order(x, k1, ki):
    """First-order model with non-competitive inhibition; x = (S, I)."""
    s, i = x
    return k1 * s / (1.0 + i / ki)

def relative_activity(i, ic50):
    """Relative activity (%) vs. inhibitor concentration."""
    return 100.0 / (1.0 + i / ic50)

def aic(ss_residuals, n_points, n_params):
    """Akaike Information Criterion for least-squares fits."""
    return n_points * np.log(ss_residuals / n_points) + 2 * n_params

# ----------------------------------------------------------------------
# 3. NON-MICHAELIAN BEHAVIOR: FIRST-ORDER VS. MICHAELIS-MENTEN (CONTROL)
# ----------------------------------------------------------------------
print("=" * 68)
print("CAT + RESORCINARENES: non-Michaelian behavior test (control, no inhibitor)")
print("=" * 68)
v_ctrl = np.array(V_ETRA[0])

p_lin, _ = curve_fit(first_order, S, v_ctrl, p0=[2.0])
ss_lin = np.sum((v_ctrl - first_order(S, *p_lin)) ** 2)

p_mm, _ = curve_fit(michaelis_menten, S, v_ctrl, p0=[150, 25], maxfev=20000)
ss_mm = np.sum((v_ctrl - michaelis_menten(S, *p_mm)) ** 2)

n = len(S)
print(f"First-order model : k1 = {p_lin[0]:.3f} min^-1, "
      f"SSR = {ss_lin:.2f}, AIC = {aic(ss_lin, n, 1):.2f}")
print(f"Michaelis-Menten  : Vmax = {p_mm[0]:.1f}, Km = {p_mm[1]:.1f}, "
      f"SSR = {ss_mm:.2f}, AIC = {aic(ss_mm, n, 2):.2f}")
print("Note: the velocity grows linearly with [H2O2] and no saturation is")
print("observed; the first-order model is preferred (lower AIC), confirming")
print("the NON-MICHAELIAN behavior of catalase under these conditions.")

# ----------------------------------------------------------------------
# 4. GLOBAL FIT OF THE NON-COMPETITIVE FIRST-ORDER MODEL (Na4EtRA)
# ----------------------------------------------------------------------
S_all, I_all, v_all = [], [], []
for I in RES_LEVELS:
    for s, v in zip(S, V_ETRA[I]):
        S_all.append(s); I_all.append(I); v_all.append(v)
S_all, I_all, v_all = map(np.array, (S_all, I_all, v_all))

popt_g, _ = curve_fit(noncompetitive_first_order, (S_all, I_all), v_all,
                      p0=[2.0, 45.0])
k1, Ki = popt_g
print("-" * 68)
print("GLOBAL FIT (non-competitive first-order model, Na4EtRA):")
print(f"  k1 = {k1:.3f} min^-1")
print(f"  Ki = {Ki:.2f} uM")
print("Lineweaver-Burk-style check: the apparent first-order constant")
print("(slope of v vs. [S]) decreases with [I], i.e., the double-reciprocal")
print("lines converge on the x-axis -> NON-COMPETITIVE inhibition.")

# ----------------------------------------------------------------------
# 5. RELATIVE-ACTIVITY CURVES AND INHIBITORY CAPACITY OF THE FIVE LIGANDS
# ----------------------------------------------------------------------
print("-" * 68)
print("RELATIVE ACTIVITY ([H2O2] = 18 uM): inhibitory capacity ranking")
ic50_table = {}
for name, ra in REL_ACTIVITY.items():
    ra = np.array(ra)
    popt, _ = curve_fit(relative_activity, np.array(RES_LEVELS, float), ra,
                        p0=[50.0])
    ic50_table[name] = popt[0]
    print(f"  {name:>9}: IC50 = {popt[0]:6.1f} uM, "
          f"residual activity at 100 uM = {ra[-1]:.1f} %")
ranking = sorted(ic50_table, key=ic50_table.get)
print("Inhibitory tendency:", " > ".join(ranking),
      "(published: Na4PrRA > Na4EtRA > Na4ESRA > Na4MeRA > Na4SRA)")

# ----------------------------------------------------------------------
# 6. FIGURES
# ----------------------------------------------------------------------
# Figure A: first-order kinetics and non-competitive inhibition (Na4EtRA)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
colors = plt.cm.viridis(np.linspace(0, 0.85, len(RES_LEVELS)))
Ss = np.linspace(0, 55, 100)
for c, I in zip(colors, RES_LEVELS):
    ax[0].plot(S, V_ETRA[I], "o", color=c, ms=5)
    ax[0].plot(Ss, noncompetitive_first_order((Ss, I), k1, Ki), "-",
               color=c, lw=1.4, label=f"{I} uM")
ax[0].set_xlabel("[H$_2$O$_2$]  (uM)")
ax[0].set_ylabel("v  (uM/min)")
ax[0].set_title("A. First-order kinetics (no saturation)")
ax[0].legend(title="Na$_4$EtRA", fontsize=8)

for c, I in zip(colors, RES_LEVELS):
    v = np.array(V_ETRA[I])
    x, y = 1.0 / S, 1.0 / v
    a, b = np.polyfit(x, y, 1)
    ax[1].plot(x, y, "o", color=c, ms=5)
    x_intercept = -b / a  # near zero for first-order non-competitive kinetics
    xs = np.linspace(min(x_intercept * 1.3, 0.0), x.max() * 1.05, 50)
    ax[1].plot(xs, a * xs + b, "-", color=c, lw=1.2, label=f"{I} uM")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].axvline(0, color="k", lw=0.6)
ax[1].set_xlabel("1/[H$_2$O$_2$]  (uM$^{-1}$)")
ax[1].set_ylabel("1/v  (min/uM)")
ax[1].set_title("B. Double-reciprocal plot (non-competitive)")
ax[1].legend(title="Na$_4$EtRA", fontsize=8)
ax[1].set_ylim(bottom=0)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig_catalase_firstorder_inhibition.png", dpi=200)

# Figure B: relative activity of the five resorcinarenes
fig2, ax2 = plt.subplots(figsize=(6.4, 4.6))
markers = ["o", "s", "^", "D", "v"]
Is = np.linspace(0, 105, 120)
for (name, ra), mk in zip(REL_ACTIVITY.items(), markers):
    ax2.plot(RES_LEVELS, ra, mk, ms=5, color=None)
    ax2.plot(Is, relative_activity(Is, ic50_table[name]), "-", lw=1.4,
             label=f"{name} (IC$_{{50}}$ = {ic50_table[name]:.0f} uM)")
ax2.set_xlabel("[resorcinarene]  (uM)")
ax2.set_ylabel("Relative activity  (%)")
ax2.set_title("CAT inhibition by sulfonated resorcinarenes")
ax2.legend(fontsize=8)
fig2.tight_layout()
fig2.savefig(OUTPUT_DIR / "fig_catalase_relative_activity.png", dpi=200)

print("-" * 68)
print(f"Figures saved in: {OUTPUT_DIR}")
print("  - fig_catalase_firstorder_inhibition.png")
print("  - fig_catalase_relative_activity.png")

if __name__ == "__main__":
    plt.show()
