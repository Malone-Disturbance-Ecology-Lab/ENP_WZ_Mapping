"""
Graphical abstract — Zhang & Malone, SESE land-cover change (1984–2024)
Repo: Malone-Disturbance-Ecology-Lab/ENP_WZ_Mapping

Data sources (all final R2 manuscript values):
  * Panel 1: annual DZ/TZ/VZ area for the whole SESE, digitized from R2 Fig. 5b
    (data/fig5b_digitized_SESE_area_km2.csv). Check: 1984->2024 change
    DZ +22.2, TZ -25.5, VZ +3.1 km^2 vs. reported +22.3 / -25.6 / +3.3; annual
    sums = 182.6-183.0 km^2 (study area 182.9 km^2).
  * Panel 2 map: SESE outline + 100 transects from the repo workspace
    (aoi, ENP.transects). Subregions by transect ID (>61 W, 37-61 C, <37 E).
  * Panel 2 bars: Table 1 mean transect rates (m/yr); subregional mean
    salinity 1991-2024 (Sec. 3.1); VZ breakpoint 25.87 PSU (Table S6).
  * Panel 3: 2074 scenario areas (Sec. 3.4 + Fig. 8 endpoints, digitized).
Swap in the final EDI raster series when available: only DATA dicts change.
"""
import pandas as pd, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly, FancyArrowPatch, FancyBboxPatch
from matplotlib.gridspec import GridSpec

mpl.rcParams.update({
    "font.family": "Fira Sans", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#8a8a8a", "axes.linewidth": 0.7,
    "xtick.color": "#555", "ytick.color": "#555",
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

# ---- colours: manuscript DZ/TZ/VZ palette (TZ deepened for contrast) ----
DZ, TZ, VZ = "#B08C60", "#7FC46E", "#00461E"
INK, INK2, MUTED = "#1d1d1b", "#4a4a48", "#8a8a88"
SAL = {"West": "#9ECAE1", "Central": "#4292C6", "East": "#08306B"}  # sequential: low -> high salinity
FW = "#2B7BB9"

# ---- data ----
ts = pd.read_csv("data/fig5b_digitized_SESE_area_km2.csv")
sub = pd.DataFrame({
    "region":   ["West", "Central", "East"],
    "salinity": [17.6, 21.9, 29.8],          # mean annual PSU, 1991-2024
    "DZ_rate":  [16.5, 32.3, 23.5],          # m/yr, Table 1
    "VZ_rate":  [19.5, -13.8, -2.5],         # m/yr, Table 1
})
scen = pd.DataFrame({   # km^2 in 2074
    "scenario": ["More freshwater\n(low-salinity rates)", "Current\ntrajectory", "Less freshwater\n(high-salinity rates)"],
    "VZ": [97.2, 76.2, 63.3],
    "DZ": [53.6, 66.6, 74.0],
})
now = {"VZ": float(ts.loc[ts.year == 2024, "VZ"].iloc[0]), "DZ": float(ts.loc[ts.year == 2024, "DZ"].iloc[0])}

aoi = pd.read_csv("data/geom_aoi.csv")
tg = pd.read_csv("data/geom_ENP.transects.csv")
ta = pd.read_csv("data/attr_ENP.transects.csv")
ta["feat"] = np.arange(1, len(ta) + 1)
ta["region"] = np.where(ta.ID > 61, "West", np.where(ta.ID >= 37, "Central", "East"))

# ---- canvas: Elsevier landscape GA (2.5:1) ----
W, H = 13.33, 5.6
fig = plt.figure(figsize=(W, H), dpi=300)
fig.patch.set_facecolor("white")

fig.text(0.022, 0.945, "Freshwater sets the balance between wetland degradation and vegetation expansion",
         fontsize=16.5, weight="bold", color=INK, va="center")
fig.text(0.022, 0.895, "Southeastern Saline Everglades (182.9 km²)  ·  41 years of Landsat, machine-learning classification "
         "with temporal transfer learning  ·  100 transects × 3 salinity subregions",
         fontsize=9.6, color=INK2, va="center")

gs = GridSpec(2, 3, figure=fig, left=0.045, right=0.985, top=0.775, bottom=0.185,
              width_ratios=[1.0, 1.25, 0.95], height_ratios=[1.05, 1], wspace=0.36, hspace=0.32)

HEAD_Y, SUB_Y = 0.848, 0.812
def step_label(ax, n, title, sub=None, **_):
    x0 = ax.get_position().x0
    fig.text(x0 - 0.004, HEAD_Y, f"{n}", fontsize=10, weight="bold", color="white", ha="center", va="center",
             bbox=dict(boxstyle="circle,pad=0.28", fc=INK, ec="none"))
    fig.text(x0 + 0.012, HEAD_Y, title, fontsize=11.5, weight="bold", color=INK, va="center")
    if sub:
        fig.text(x0 + 0.012, SUB_Y, sub, fontsize=8.4, color=INK2, va="center")

# ================= Panel 1: 41-yr time series =================
a1 = fig.add_subplot(gs[:, 0])
a1.axvspan(2014, 2024.6, color=FW, alpha=0.08, lw=0)
a1.text(2019.3, 113.5, "↑ freshwater flow\n(restoration, since 2014)", ha="center", va="top",
        fontsize=7.6, color=FW, linespacing=1.05)
for col, c, lw in [("TZ", TZ, 1.6), ("VZ", VZ, 1.6), ("DZ", DZ, 2.2)]:
    a1.plot(ts.year, ts[col], color=c, lw=lw, solid_capstyle="round", zorder=3)
    a1.scatter(ts.year, ts[col], s=7, color=c, edgecolor="white", linewidth=0.5, zorder=4)
# VZ trend since 2014 (paper: 0.3 km2/yr)
seg = ts[ts.year >= 2014]
b = np.polyfit(seg.year, seg.VZ, 1)
xx = np.array([2014, 2024]); a1.plot(xx, np.polyval(b, xx), color="#9a9a9a", lw=1.6, zorder=2, ls=(0, (3, 1.5)))

lab = dict(fontsize=8.8, va="center", ha="left", clip_on=False)
a1.text(2025.6, 69.9, "Transition", color=INK, weight="bold", **lab)
a1.text(2025.6, 64.5, "−25.6 km²", color=INK2, **lab)
a1.text(2025.6, 57.0, "Vegetated", color=INK, weight="bold", **lab)
a1.text(2025.6, 52.1, "+3.3 km²", color=INK2, **lab)
a1.text(2025.6, 44.0, "Degraded", color=INK, weight="bold", **lab)
a1.text(2025.6, 39.0, "+22.3 km²", color=DZ, weight="bold", **lab)
for y, c in [(67.2, TZ), (54.6, VZ), (41.5, DZ)]:
    a1.plot([2025.0], [y], marker="s", ms=5, color=c, clip_on=False)
a1.set_xlim(1983, 2024.8); a1.set_ylim(18, 115)
a1.set_xticks([1984, 1994, 2004, 2014, 2024]); a1.set_yticks([20, 40, 60, 80, 100])
a1.set_ylabel("Area (km²)", color=INK2, fontsize=9)
a1.grid(axis="y", color="#ececec", lw=0.6); a1.set_axisbelow(True)
a1.text(2004, 18.8, "Degraded zone (DZ) expands landward 0.8 km", fontsize=7.5, color=DZ, ha="center", va="bottom")
step_label(a1, 1, "Degradation spreads", "Annual land-cover area, 1984–2024", y=1.075)

# ================= Panel 2a: map =================
a2 = fig.add_subplot(gs[0, 1])
x, y = aoi.x.values, aoi.y.values
for f, g in tg.groupby("feat"):
    r = ta.loc[ta.feat == f, "region"].iloc[0]
    a2.plot(g.x, g.y, color=SAL[r], lw=1.15, solid_capstyle="round", zorder=2)
a2.set_aspect(1 / np.cos(np.deg2rad(25.25)))
pad = 0.012
a2.set_xlim(tg.x.min() - 0.05, tg.x.max() + pad); a2.set_ylim(tg.y.min() - 0.045, tg.y.max() + 0.03)
a2.axis("off")
# subregion labels (salinity)
cent = ta.merge(tg.groupby("feat")[["x", "y"]].mean().reset_index(), on="feat").groupby("region")[["x", "y"]].mean()
offs = {"West": (0.0, -0.062), "Central": (0.01, -0.06), "East": (0.035, -0.055)}
for r in ["West", "Central", "East"]:
    s = sub.set_index("region").loc[r, "salinity"]
    cx, cy = cent.loc[r]
    a2.text(cx + offs[r][0], cy + offs[r][1], f"{r}\n{s:.1f} PSU", ha="center", va="center", fontsize=8.6,
            color="white" if r != "West" else INK, weight="bold", linespacing=1.0,
            bbox=dict(boxstyle="round,pad=0.3,rounding_size=0.4", fc=SAL[r], ec="none"), zorder=5)
# freshwater input arrow (Taylor Slough, west)
wx, wy = tg[tg.feat.isin(ta.loc[ta.region == "West", "feat"])][["x", "y"]].agg(["min", "max"]).values
a2.add_patch(FancyArrowPatch((wx[0] - 0.035, wy[1] + 0.03), (wx[0] + 0.01, wy[1] - 0.005),
             arrowstyle="-|>,head_length=4,head_width=2.6", lw=1.6, color=FW, zorder=6,
             connectionstyle="arc3,rad=0.25", clip_on=False))
a2.text(wx[0] - 0.04, wy[1] + 0.035, "Freshwater input\n(Taylor Slough)", fontsize=7.6, color=FW,
        ha="left", va="bottom", linespacing=1.05, clip_on=False)
a2.text(tg.x.mean() + 0.02, tg.y.min() - 0.035, "Florida Bay", fontsize=8, color="#7d9fb8", style="italic", ha="center")
step_label(a2, 2, "Salinity steers the outcome", "100 transects across a west–east salinity gradient", y=1.20, dy=0.155)

# ================= Panel 2b: rates by subregion =================
a3 = fig.add_subplot(gs[1, 1])
xi = np.arange(3); bw = 0.34
bd = a3.bar(xi - bw / 2 - 0.02, sub.DZ_rate, bw, color=DZ, zorder=3)
bv = a3.bar(xi + bw / 2 + 0.02, sub.VZ_rate, bw, color=VZ, zorder=3)
a3.axhline(0, color="#777", lw=0.8, zorder=4)
for rect, v in list(zip(bd, sub.DZ_rate)) + list(zip(bv, sub.VZ_rate)):
    a3.text(rect.get_x() + rect.get_width() / 2, v + (1.3 if v >= 0 else -1.3), f"{v:+.1f}",
            ha="center", va="bottom" if v >= 0 else "top", fontsize=7.8, color=INK2)
a3.set_xticks(xi)
a3.set_xticklabels([f"{r}\n{s:.1f} PSU" for r, s in zip(sub.region, sub.salinity)], fontsize=8.2, color=INK2)
for t, r in zip(a3.get_xticklabels(), sub.region):
    t.set_color(SAL[r] if r != "West" else "#5d9cc6")
a3.tick_params(axis="x", length=0)
a3.set_ylim(-22, 42); a3.set_yticks([-20, 0, 20, 40])
a3.set_ylabel("Rate (m yr⁻¹)", color=INK2, fontsize=8.6)
a3.grid(axis="y", color="#ececec", lw=0.6); a3.set_axisbelow(True)
a3.legend([bd, bv], ["Degraded zone", "Vegetated zone"], frameon=False, fontsize=7.8, ncol=2,
          loc="upper left", bbox_to_anchor=(0.0, 1.17), handlelength=1.0, columnspacing=1.0)
a3.text(2.55, -10.0, "Dense vegetation declines\nabove ~26 PSU (breakpoint)", fontsize=7.6, color=VZ, ha="right", va="top",
        linespacing=1.05)
a3.text(-0.48, 39, "expansion ↑", fontsize=7, color=MUTED, va="center")
a3.text(-0.48, -17.5, "retraction ↓", fontsize=7, color=MUTED, va="center")
a3.set_xlim(-0.55, 2.55)

# ================= Panel 3: 2074 scenarios =================
a4 = fig.add_subplot(gs[:, 2])
yi = np.arange(3)[::-1] * 1.0
bh = 0.32
for i, row in scen.iterrows():
    yv, yd = yi[i] + bh / 2 + 0.02, yi[i] - bh / 2 - 0.02
    a4.barh(yv, row.VZ, bh, color=VZ, zorder=3)
    a4.barh(yd, row.DZ, bh, color=DZ, zorder=3)
    a4.text(row.VZ + 1.5, yv, f"{row.VZ:.0f}", va="center", fontsize=8.4, color=INK, weight="bold")
    a4.text(row.DZ + 1.5, yd, f"{row.DZ:.0f}", va="center", fontsize=8.4, color=INK2)
for k, c in [("VZ", VZ), ("DZ", DZ)]:
    a4.axvline(now[k], color=c, lw=1.0, ls=(0, (2, 2)), zorder=4, ymin=0.02, ymax=0.9)
a4.text(now["VZ"] - 0.8, 2.72, f"2024\nVZ", fontsize=7, color=VZ, ha="right", va="bottom", linespacing=1.0)
a4.text(now["DZ"] - 0.8, 2.72, f"2024\nDZ", fontsize=7, color=DZ, ha="right", va="bottom", linespacing=1.0)
a4.set_yticks(yi); a4.set_yticklabels(scen.scenario, fontsize=8.6, color=INK, linespacing=1.05)
for t, c in zip(a4.get_yticklabels(), [FW, INK2, "#a0522d"]):
    t.set_color(c)
a4.tick_params(axis="y", length=0)
a4.set_xlim(0, 112); a4.set_ylim(-1.0, 3.05); a4.set_xticks([0, 25, 50, 75, 100])
a4.set_xlabel("Area in 2074 (km²)", color=INK2, fontsize=8.6)
a4.grid(axis="x", color="#ececec", lw=0.6); a4.set_axisbelow(True)
a4.spines["left"].set_visible(False)
a4.legend(handles=[mpl.patches.Patch(color=VZ, label="Vegetated zone"), mpl.patches.Patch(color=DZ, label="Degraded zone")],
          frameon=False, fontsize=7.8, ncol=2, loc="lower right", bbox_to_anchor=(1.0, 0.0), handlelength=1.0)
step_label(a4, 3, "Freshwater shifts the future", "50-year land-cover scenarios", y=1.075)

# ================= Takeaway band =================
fig.patches.append(FancyBboxPatch((0.022, 0.02), 0.963, 0.066, boxstyle="round,pad=0,rounding_size=0.008",
                                  transform=fig.transFigure, fc="#eaf2f8", ec="none", zorder=0))
fig.text(0.035, 0.053, "Takeaway", fontsize=10, weight="bold", color=FW, va="center")
fig.text(0.105, 0.053, "Where freshwater input is highest, degradation is slowest and dense vegetation expands — "
         "increasing freshwater flow is a lever for resilient coastal wetlands under sea-level rise.",
         fontsize=10, color=INK, va="center")

for ext in ["png", "pdf", "svg"]:
    fig.savefig(f"out/Graphical_Abstract_SESE.{ext}", dpi=300, facecolor="white")
fig.savefig("out/Graphical_Abstract_SESE.tiff", dpi=300, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
print("saved")
