"""
Umerkot / Nara Canal — agricultural land suitability & its trend (2000–2025)
STEP 2: model + trend analysis on the layers exported by gee/01_extract_data.js

Usage:  python run_analysis.py <folder with the .npy files> <output folder>

What it does
  1. Loads yearly layers (peak NDVI, dry-season salinity & brightness indices,
     rainfall, land-surface temperature) and static layers (slope, soil sand %,
     land cover, surface water).
  2. Turns each into a 0–1 "goodness" score (fuzzy linear membership).
  3. Weights them with AHP (Analytic Hierarchy Process) -> yearly suitability 0–1.
  4. Classifies current suitability (mean 2023–25) into FAO classes S1/S2/S3/N.
  5. Per pixel: Mann-Kendall test + Sen's slope of suitability 2000–2025,
     plus a rainfall-adjusted regression, to find where suitability is
     falling fastest and whether that decline is real (not just drought).
  6. Writes maps (PNG), GIS rasters (ESRI ASCII .asc), and summary tables (CSV).
"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, TwoSlopeNorm
from scipy import stats

# ------------------------------------------------------------------ grid
WEST, NORTH, PIX = 69.3, 25.6, 0.001
NX, NY = 800, 700
EXTENT = [WEST, WEST + NX * PIX, NORTH - NY * PIX, NORTH]
YEARS = np.arange(2000, 2026)
NODATA = -32768

TOWNS = {  # approximate town locations, for orientation only
    "Umerkot": (69.736, 25.362), "Kunri": (69.566, 25.178),
    "Samaro": (69.394, 25.282), "Pithoro": (69.380, 25.510),
    "Chhor": (69.780, 25.510),
}

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
fig_dir = os.path.join(out, "figures"); os.makedirs(fig_dir, exist_ok=True)
gis_dir = os.path.join(out, "gis"); os.makedirs(gis_dir, exist_ok=True)


def load_stack(name, prefix, scale):
    """Load an exported NPY (structured array) into float array [years, y, x]."""
    a = np.load(os.path.join(src, name + ".npy"))
    layers = []
    for y in YEARS:
        b = a[f"{prefix}_{y}"].astype("float32")
        b[b == NODATA] = np.nan
        layers.append(b / scale)
    return np.stack(layers)


ndvi = load_stack("ndvi_p90", "ndvi_p90", 10000)
ndvi_dry = load_stack("ndvi_dry", "ndvi_dry", 10000)
si = load_stack("si_dry", "si_dry", 10000)
bi = load_stack("bi_dry", "bi_dry", 10000)
rain = load_stack("rain_mm", "rain", 1)
lst = load_stack("lst_c_x100", "lst", 100)
nsc = load_stack("n_scenes", "n", 1)

st = np.load(os.path.join(src, "static.npy"))
def sfield(f):
    b = st[f].astype("float32"); b[b == NODATA] = np.nan; return b
slope, sand, clay = sfield("slope_deg"), sfield("sand_pct"), sfield("clay_pct")
wc, water, elev = sfield("worldcover"), sfield("water_occ"), sfield("elev_m")
rgb = np.load(os.path.join(src, "rgb_2025.npy"))
rgb = np.dstack([rgb["red"], rgb["green"], rgb["blue"]])

# ------------------------------------------------------------------ masks & zones
# Exclude permanent water and built-up land from the analysis
exclude = (wc == 80) | (wc == 50) | (water > 50)
# "Farmed zone": ESA cropland OR clearly green (peak NDVI >= 0.4) in at least 3 years
ever_green = (np.nansum(ndvi >= 0.40, axis=0) >= 3)
farmed = ((wc == 40) | ever_green) & ~exclude
desert = ~farmed & ~exclude
valid = ~exclude


def ramp(x, lo, hi, increasing=True):
    """Fuzzy linear membership: 0 at 'bad' end, 1 at 'good' end."""
    z = np.clip((x - lo) / (hi - lo), 0, 1)
    return z if increasing else 1 - z


# ------------------------------------------------------------------ factors
# Salinity proxy (0 = no salt signal, 1 = strong): mean of rescaled SI and BI,
# computed only where the dry-season surface is bare (NDVI_dry < 0.25);
# vegetated dry-season pixels are taken to carry little surface salt.
si_n = ramp(si, np.nanpercentile(si, 5), np.nanpercentile(si, 95))
bi_n = ramp(bi, np.nanpercentile(bi, 5), np.nanpercentile(bi, 95))
bare = ramp(ndvi_dry, 0.35, 0.15)          # 1 when bare, 0 when well-vegetated
salinity = (0.5 * si_n + 0.5 * bi_n) * bare

F = {
    "Low salinity":      1 - salinity,
    "Crop vigour (NDVI)": ramp(ndvi, 0.10, 0.70),
    "Soil texture":      np.broadcast_to(ramp(sand, 70, 30), ndvi.shape),   # less sand = better (local range 25–70 %)
    "Rainfall":          ramp(rain, 100, 500),
    "Cooler surface":    ramp(lst, 50, 38),                                 # hotter = worse
    "Gentle slope":      np.broadcast_to(ramp(slope, 8, 1), ndvi.shape),
}

# ------------------------------------------------------------------ AHP weights
names = list(F)
# Pairwise importance (Saaty 1–9 scale), row factor vs column factor
P = np.array([
    # Sal  NDVI Soil Rain Heat Slope
    [1,    2,   3,   4,   4,   6],   # Low salinity
    [1/2,  1,   2,   3,   3,   5],   # Crop vigour
    [1/3,  1/2, 1,   2,   2,   4],   # Soil texture
    [1/4,  1/3, 1/2, 1,   1,   3],   # Rainfall
    [1/4,  1/3, 1/2, 1,   1,   3],   # Cooler surface
    [1/6,  1/5, 1/4, 1/3, 1/3, 1],   # Gentle slope
])
eigval, eigvec = np.linalg.eig(P)
k = np.argmax(eigval.real)
w = np.abs(eigvec[:, k].real); w /= w.sum()
lam = eigval[k].real; n = len(P)
CI = (lam - n) / (n - 1); RI = 1.24   # Saaty random index, n=6
CR = CI / RI
weights = dict(zip(names, w))

# ------------------------------------------------------------------ suitability
suit = sum(weights[k] * F[k] for k in names)       # [years, y, x]
suit[:, ~valid] = np.nan
current = np.nanmean(suit[-3:], axis=0)            # 2023–2025 average
early = np.nanmean(suit[:3], axis=0)               # 2000–2002 average

# FAO-style classes
bins = [0.70, 0.55, 0.40]
cls = np.full(current.shape, np.nan)
cls[current >= 0.70] = 1                                # S1 highly suitable
cls[(current >= 0.55) & (current < 0.70)] = 2           # S2 moderately suitable
cls[(current >= 0.40) & (current < 0.55)] = 3           # S3 marginally suitable
cls[current < 0.40] = 4                                 # N  not suitable
cls[~valid] = np.nan

# ------------------------------------------------------------------ trend tests
def mann_kendall_sen(Y, t):
    """Vectorised Mann-Kendall (no-ties variance) + Sen's slope.
    Y: [T, N] with NaNs allowed. Returns tau-like S/z, p, sen slope."""
    T, N = Y.shape
    i, j = np.triu_indices(T, 1)
    d = Y[j] - Y[i]                                  # [pairs, N]
    dt = (t[j] - t[i])[:, None].astype("float32")
    ok = ~np.isnan(d)
    S = np.nansum(np.sign(d), axis=0)
    n_obs = np.sum(~np.isnan(Y), axis=0).astype("float64")
    var = n_obs * (n_obs - 1) * (2 * n_obs + 5) / 18.0
    z = np.where(S > 0, (S - 1) / np.sqrt(var), np.where(S < 0, (S + 1) / np.sqrt(var), 0))
    p = 2 * (1 - stats.norm.cdf(np.abs(z)))
    slopes = np.where(ok, d / dt, np.nan)
    sen = np.nanmedian(slopes, axis=0)
    p[n_obs < 10] = np.nan; sen[n_obs < 10] = np.nan
    return z, p, sen


def ols_trend_with_rain(Y, t, R):
    """Per-pixel OLS: Y = a + b*year + c*rain.  Returns b and its p-value."""
    T, N = Y.shape
    ok = ~np.isnan(Y) & ~np.isnan(R)
    tc = (t - t.mean()).astype("float64")[:, None] * np.ones((1, N))
    X = np.stack([np.ones_like(tc), tc, R.astype("float64")], axis=-1)  # [T,N,3]
    X = np.where(ok[..., None], X, 0.0); Yz = np.where(ok, Y, 0.0)
    XtX = np.einsum("tni,tnj->nij", X, X)
    XtY = np.einsum("tni,tn->ni", X, Yz)
    nobs = ok.sum(0)
    good = (nobs >= 10) & (np.abs(np.linalg.det(XtX)) > 1e-9)
    b = np.full(N, np.nan); pval = np.full(N, np.nan)
    XtXg = XtX[good]
    beta = np.linalg.solve(XtXg, XtY[good][..., None])[..., 0]
    resid = Yz[:, good] - np.einsum("tni,ni->tn", X[:, good], beta)
    resid = np.where(ok[:, good], resid, 0)
    dof = nobs[good] - 3
    s2 = (resid ** 2).sum(0) / dof
    inv = np.linalg.inv(XtXg)
    se = np.sqrt(s2 * inv[:, 1, 1])
    tstat = beta[:, 1] / se
    b[good] = beta[:, 1]
    pval[good] = 2 * (1 - stats.t.cdf(np.abs(tstat), dof))
    return b, pval


flat = lambda a: a.reshape(a.shape[0], -1)
idx = np.where(valid.ravel())[0]

z = np.full(NX * NY, np.nan); p = z.copy(); sen = z.copy()
b_rain = z.copy(); p_rain = z.copy()
sal_sen = z.copy(); sal_p = z.copy()
for chunk in np.array_split(idx, 40):
    Ys = flat(suit)[:, chunk]
    zz, pp, ss = mann_kendall_sen(Ys, YEARS)
    z[chunk], p[chunk], sen[chunk] = zz, pp, ss
    bb, pb = ols_trend_with_rain(Ys, YEARS, flat(rain)[:, chunk])
    b_rain[chunk], p_rain[chunk] = bb, pb
    _, sp, sz = mann_kendall_sen(flat(salinity)[:, chunk], YEARS)
    sal_sen[chunk], sal_p[chunk] = sz, sp
z, p, sen, b_rain, p_rain, sal_sen, sal_p = [a.reshape(NY, NX) for a in (z, p, sen, b_rain, p_rain, sal_sen, sal_p)]

# Decline categories (per decade change in the 0–1 score)
sen_dec = sen * 10
trend_cls = np.full(sen.shape, np.nan)
sig = p < 0.05
robust = sig & (p_rain < 0.05) & (b_rain < 0)
trend_cls[valid] = 3                                     # no significant change
trend_cls[sig & (sen_dec > 0)] = 4                       # improving
trend_cls[sig & (sen_dec < 0)] = 2                       # declining
trend_cls[robust & (sen_dec <= -0.04)] = 1               # fast decline, holds after rain adjustment

# ------------------------------------------------------------------ robustness checks
# (1) Landsat 9 (from late 2021) roughly doubles the number of clear scenes per year,
#     which can nudge "peak NDVI" upward. Re-run the trend on 2000–2021 only.
# (2) Is "improving" still improving once rainfall is accounted for?
sub = YEARS <= 2021
sen21 = np.full(NX * NY, np.nan); p21 = sen21.copy()
for chunk in np.array_split(idx, 40):
    _, pp, ss = mann_kendall_sen(flat(suit[sub])[:, chunk], YEARS[sub])
    sen21[chunk], p21[chunk] = ss, pp
sen21, p21 = sen21.reshape(NY, NX), p21.reshape(NY, NX)
fast = trend_cls == 1
impr = trend_cls == 4
robust_info = {
    "fast_decline_pixels_also_negative_2000_2021_pct": round(100 * float(np.nanmean(sen21[fast] < 0)), 1),
    "fast_decline_pixels_also_significant_2000_2021_pct": round(100 * float(np.nanmean((sen21[fast] < 0) & (p21[fast] < 0.05))), 1),
    "improving_pixels_still_improving_after_rain_adjustment_pct": round(100 * float(np.nanmean((b_rain[impr] > 0) & (p_rain[impr] < 0.05))), 1),
    "improving_pixels_also_positive_2000_2021_pct": round(100 * float(np.nanmean(sen21[impr] > 0)), 1),
    "farmed_zone_mean_trend_per_decade_2000_2025": round(float(np.nanmean(sen_dec[farmed])), 4),
    "farmed_zone_mean_trend_per_decade_2000_2021": round(float(np.nanmean(sen21[farmed] * 10)), 4),
}

# ------------------------------------------------------------------ write GIS rasters
def write_asc(arr, name, nodata=-9999):
    a = np.where(np.isnan(arr), nodata, arr)
    hdr = (f"ncols {NX}\nnrows {NY}\nxllcorner {WEST}\nyllcorner {NORTH - NY*PIX:.6f}\n"
           f"cellsize {PIX}\nNODATA_value {nodata}\n")
    with open(os.path.join(gis_dir, name + ".asc"), "w") as f:
        f.write(hdr); np.savetxt(f, a, fmt="%.5g")
    with open(os.path.join(gis_dir, name + ".prj"), "w") as f:
        f.write('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],'
                'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]')

write_asc(current, "suitability_2023_25")
write_asc(cls, "suitability_class_2023_25")
write_asc(sen_dec, "suitability_trend_per_decade")
write_asc(p, "suitability_trend_pvalue")
write_asc(trend_cls, "trend_class")
write_asc(np.nanmean(salinity[-3:], 0), "salinity_proxy_2023_25")
write_asc(sal_sen * 10, "salinity_trend_per_decade")

# ------------------------------------------------------------------ summaries
pix_km2 = (PIX * 111.32 * np.cos(np.radians(25.25))) * (PIX * 110.6)
rows = []
for zone, m in [("All land", valid), ("Farmed zone", farmed), ("Desert / rangeland", desert)]:
    r = {"zone": zone, "area_km2": round(m.sum() * pix_km2)}
    for c, lab in [(1, "S1"), (2, "S2"), (3, "S3"), (4, "N")]:
        r[f"{lab}_pct"] = round(100 * np.nansum((cls == c) & m) / m.sum(), 1)
    for c, lab in [(1, "fast_decline"), (2, "declining"), (3, "no_change"), (4, "improving")]:
        r[f"{lab}_pct"] = round(100 * np.nansum((trend_cls == c) & m) / m.sum(), 1)
    r["mean_suit_2000_02"] = round(float(np.nanmean(early[m])), 3)
    r["mean_suit_2023_25"] = round(float(np.nanmean(current[m])), 3)
    rows.append(r)
import csv
with open(os.path.join(out, "summary_by_zone.csv"), "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)

# Yearly regional means
ts = []
for i, y in enumerate(YEARS):
    ts.append({"year": int(y),
               "suitability_farmed": float(np.nanmean(suit[i][farmed])),
               "suitability_desert": float(np.nanmean(suit[i][desert])),
               "salinity_farmed": float(np.nanmean(salinity[i][farmed])),
               "peak_ndvi_farmed": float(np.nanmean(ndvi[i][farmed])),
               "rain_mm": float(np.nanmean(rain[i][valid])),
               "lst_c": float(np.nanmean(lst[i][valid])),
               "median_clear_scenes": float(np.nanmedian(nsc[i]))})
with open(os.path.join(out, "yearly_means.csv"), "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(ts[0])); wr.writeheader(); wr.writerows(ts)

# Hotspot table: 5 x 5 km blocks ranked by share of fast-declining land
blk = 50
hot = []
for by in range(0, NY, blk):
    for bx in range(0, NX, blk):
        tc = trend_cls[by:by+blk, bx:bx+blk]; v = valid[by:by+blk, bx:bx+blk]
        if v.sum() < 0.5 * blk * blk: continue
        lon = WEST + (bx + blk / 2) * PIX; lat = NORTH - (by + blk / 2) * PIX
        near = min(TOWNS, key=lambda t: (TOWNS[t][0]-lon)**2 + (TOWNS[t][1]-lat)**2)
        dkm = np.hypot((TOWNS[near][0]-lon)*100.7, (TOWNS[near][1]-lat)*110.6)
        hot.append({"block_centre_lat": round(lat, 3), "block_centre_lon": round(lon, 3),
                    "nearest_town": near, "km_from_town": round(dkm, 1),
                    "fast_decline_pct": round(100*np.nansum(tc == 1)/v.sum(), 1),
                    "any_decline_pct": round(100*np.nansum(tc <= 2)/v.sum(), 1),
                    "mean_trend_per_decade": round(float(np.nanmean(sen_dec[by:by+blk, bx:bx+blk])), 4),
                    "farmed_pct": round(100*farmed[by:by+blk, bx:bx+blk].sum()/v.sum(), 1)})
hot.sort(key=lambda r: -r["fast_decline_pct"])
with open(os.path.join(out, "decline_hotspots_5km.csv"), "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(hot[0])); wr.writeheader(); wr.writerows(hot)

meta = {"ahp_weights": {k: round(float(v), 3) for k, v in weights.items()},
        "ahp_consistency_ratio": round(float(CR), 3),
        "pixel_km2": round(pix_km2, 4),
        "farmed_km2": round(float(farmed.sum() * pix_km2)),
        "desert_km2": round(float(desert.sum() * pix_km2)),
        "class_thresholds": {"S1": ">=0.70", "S2": "0.55-0.70", "S3": "0.40-0.55", "N": "<0.40"},
        "fast_decline_rule": "MK p<0.05 AND rain-adjusted OLS slope<0 with p<0.05 AND Sen slope <= -0.04 per decade",
        "robustness": robust_info}
json.dump(meta, open(os.path.join(out, "model_info.json"), "w"), indent=2)

# ------------------------------------------------------------------ figures
def base(ax, title):
    ax.set_title(title, fontsize=11, loc="left")
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")
    for t, (lo, la) in TOWNS.items():
        ax.plot(lo, la, "o", ms=4, mfc="white", mec="black")
        ax.text(lo + 0.01, la + 0.01, t, fontsize=8,
                bbox=dict(fc="white", ec="none", alpha=0.7, pad=1))

def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, name), dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 7))
ax.imshow(rgb, extent=EXTENT); base(ax, "Study area — true colour, dry season 2024–25 (Landsat)")
save(fig, "01_true_colour.png")

fig, ax = plt.subplots(figsize=(8, 7))
zone_img = np.full(current.shape, np.nan); zone_img[farmed] = 1; zone_img[desert] = 2; zone_img[exclude] = 3
ax.imshow(zone_img, extent=EXTENT, cmap=ListedColormap(["#4c9a2a", "#e3c88f", "#6fa8dc"]), vmin=0.5, vmax=3.5, interpolation="nearest")
base(ax, "Zones used in the analysis")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#4c9a2a", label="Farmed zone"), Patch(color="#e3c88f", label="Desert / rangeland"),
                   Patch(color="#6fa8dc", label="Water / built-up (excluded)")], loc="lower right", fontsize=8)
save(fig, "02_zones.png")

fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(current, extent=EXTENT, cmap="RdYlGn", vmin=0.2, vmax=0.85)
plt.colorbar(im, ax=ax, shrink=0.8, label="Suitability score (0–1)")
base(ax, "Agricultural suitability score, 2023–25 average"); save(fig, "03_suitability_score.png")

fig, ax = plt.subplots(figsize=(8, 7))
cm = ListedColormap(["#1a9641", "#a6d96a", "#fdae61", "#d7191c"])
ax.imshow(cls, extent=EXTENT, cmap=cm, vmin=0.5, vmax=4.5, interpolation="nearest")
ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(cm.colors,
          ["S1 Highly suitable", "S2 Moderately suitable", "S3 Marginally suitable", "N Not suitable"])],
          loc="lower right", fontsize=8)
base(ax, "Suitability classes (FAO style), 2023–25"); save(fig, "04_suitability_classes.png")

fig, ax = plt.subplots(figsize=(8, 7))
sal_now = np.nanmean(salinity[-3:], 0)
ax.imshow(np.where(farmed, np.nan, 1.0), extent=EXTENT, cmap=ListedColormap(["#d9d9d9"]), interpolation="nearest")
im = ax.imshow(np.where(farmed, sal_now, np.nan), extent=EXTENT, cmap="magma_r", vmin=0, vmax=0.8, interpolation="nearest")
plt.colorbar(im, ax=ax, shrink=0.8, label="Salinity proxy (0 = none, 1 = strong)")
ax.legend(handles=[Patch(color="#d9d9d9", label="Desert sand — index not reliable here")], loc="lower right", fontsize=8)
base(ax, "Surface salinity signal in the farmed zone, dry season 2023–25"); save(fig, "05_salinity_proxy.png")

fig, ax = plt.subplots(figsize=(8, 7))
v = np.where(sig, sen_dec, np.nan)
im = ax.imshow(v, extent=EXTENT, cmap="RdBu", norm=TwoSlopeNorm(0, -0.12, 0.12))
plt.colorbar(im, ax=ax, shrink=0.8, label="Change in suitability per decade (only p < 0.05 shown)")
base(ax, "Trend in suitability 2000–2025 (Sen's slope)"); save(fig, "06_trend_sen_slope.png")

fig, ax = plt.subplots(figsize=(8, 7))
cm2 = ListedColormap(["#b2182b", "#f4a582", "#f0f0f0", "#4393c3"])
ax.imshow(trend_cls, extent=EXTENT, cmap=cm2, vmin=0.5, vmax=4.5, interpolation="nearest")
ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(cm2.colors,
          ["Fast decline (robust to rainfall)", "Declining", "No significant change", "Improving"])],
          loc="lower right", fontsize=8)
for r in hot[:5]:
    ax.add_patch(plt.Rectangle((r["block_centre_lon"]-0.025, r["block_centre_lat"]-0.025), 0.05, 0.05, fill=False, ec="black", lw=1.5))
base(ax, "Where suitability is falling (boxes = top-5 hotspots)"); save(fig, "07_decline_hotspots.png")

fig, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
yy = [r["year"] for r in ts]
axs[0].plot(yy, [r["suitability_farmed"] for r in ts], "-o", ms=3, label="Farmed zone")
axs[0].plot(yy, [r["suitability_desert"] for r in ts], "-o", ms=3, label="Desert / rangeland")
axs[0].set_ylabel("Mean suitability"); axs[0].legend(fontsize=8)
axs[1].plot(yy, [r["salinity_farmed"] for r in ts], "-o", ms=3, color="#7b3294")
axs[1].set_ylabel("Salinity proxy\n(farmed zone)")
axs[2].bar(yy, [r["rain_mm"] for r in ts], color="#4393c3")
axs[2].set_ylabel("Rainfall (mm/yr)")
axs[0].set_title("Regional averages per year", loc="left")
save(fig, "08_time_series.png")

print(json.dumps(meta, indent=2))
for r in rows: print(r)
print("Top hotspots:"); [print(h) for h in hot[:8]]
