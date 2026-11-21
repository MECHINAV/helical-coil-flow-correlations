"""Coiled/straight friction ratio vs. Reynolds number for several coil curvatures."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from coil_flow import friction_coiled_ito, friction_straight

OUT = Path("results")
OUT.mkdir(exist_ok=True)

d = 0.01
re = np.linspace(5e3, 6e4, 100)
fig, ax = plt.subplots(figsize=(6, 4))
for big_d in (0.2, 0.4, 0.8):
    ratio = friction_coiled_ito(re, d, big_d) / friction_straight(re)
    ax.plot(re, ratio, label=f"D/d = {big_d / d:.0f}")
    print(f"D/d={big_d / d:>4.0f}  f_c/f_s at Re=3e4: {np.interp(3e4, re, ratio):.2f}")
ax.set_xlabel("Reynolds number")
ax.set_ylabel("f_coil / f_straight")
ax.set_title("Curvature penalty on friction factor (Ito)")
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "friction_ratio.png", dpi=150)
