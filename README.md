# Farmland suitability and salinity in Umerkot / Nara Canal, Sindh (2000–2025)

This is a remote sensing project. It uses 26 years of free satellite data to answer two questions about a farming area on the edge of the Thar Desert:

1. **How suitable is the land for farming right now?**
2. **Where is that suitability falling fastest?**

> **New to this?** Read sections 1–4 first. The glossary at the end explains every technical word.

---

## 1. The study area

- **Box:** 69.3–70.1 °E, 24.9–25.6 °N, about 80 × 78 km. That's 6,140 km² of land once water and towns are taken out.
- **What's in it:** the western half is green irrigated farmland fed by the **Nara Canal** (around Pithoro, Samaro and Kunri). The eastern half is **Thar Desert** sand dunes. **Umerkot** town sits right on the border between the two (see `outputs/figures/01_true_colour.png`).
- **Why here and not central Tharparkar?** Salinity from irrigation (salt left behind as water evaporates) is what satellites detect well, and it mostly happens on canal-irrigated land. Central Tharparkar is mostly sand, and its salinity problem is salty groundwater, which satellites can't see. This box covers both kinds of land, so the two can be compared.

## 2. How the project works (the simple version)

```
 Satellites ──► Google Earth Engine ──► 9 data files ──► Python model ──► maps + tables
 (Landsat,       (cleans clouds,         (.npy, one       (scores land,     (outputs/)
  MODIS,          makes one image        per layer)        finds trends)
  CHIRPS...)      per year)
```

**Step 1: get the data** (`gee/01_extract_data.js`, runs in Google Earth Engine)
For every year from 2000 to 2025 it builds one cloud-free image of the area from all the Landsat 5/7/8/9 photos taken that year. From those images it calculates:

| Layer | What it tells us | Source |
|---|---|---|
| Peak NDVI | How green the crops get at their best that year | Landsat |
| Dry-season Salinity Index (SI) and Brightness Index (BI) | Salt crusts are bright white on bare soil in April–June | Landsat |
| Dry-season NDVI | Whether the ground is bare or covered then | Landsat |
| Rainfall | Total rain per year | CHIRPS |
| Land surface temperature | How hot the ground gets, April–June | MODIS |
| Slope, sand %, land cover, surface water | Fixed background information | SRTM, OpenLandMap, ESA WorldCover, JRC |

**Step 2: score the land** (`analysis/run_analysis.py`)
1. Each layer is turned into a **0–1 "goodness" score**. For example, less salt means a score closer to 1.
2. The six scores are combined with **weights**, using a standard method called AHP:

   | Factor | Weight |
   |---|---|
   | Low salinity | 38 % |
   | Crop vigour (NDVI) | 25 % |
   | Soil texture (less sand) | 15 % |
   | Rainfall | 9 % |
   | Cooler surface | 9 % |
   | Gentle slope | 4 % |

   AHP consistency ratio = 0.017. Anything below 0.10 means the weights are logically consistent.
3. That gives a **suitability score from 0 to 1 for every 100 m pixel, for every year.**
4. The average for 2023–25 is sorted into FAO-style classes:
   - **S1** highly suitable: ≥ 0.70
   - **S2** moderately suitable: 0.55–0.70
   - **S3** marginally suitable: 0.40–0.55
   - **N** not suitable: < 0.40

**Step 3: find where it's getting worse**
- For each pixel, the **Mann–Kendall test** asks whether the score has *really* gone up or down over 26 years, or is just bouncing around.
- **Sen's slope** says *how fast* it's changing.
- A second test removes the effect of rainfall, so a few dry years don't get mistaken for land damage.
- A pixel is called **"fast decline"** only if it passes all of these:
  - the decline is statistically significant (p < 0.05)
  - it's still declining after the rainfall adjustment
  - it's losing at least 0.04 points per decade

## 3. Results

### Current suitability (2023–25)

| Zone | Area | S1 | S2 | S3 | N |
|---|---|---|---|---|---|
| Farmed (canal) zone | 3,289 km² | **52 %** | 22 % | 18 % | 9 % |
| Desert / rangeland | 2,851 km² | 2 % | 11 % | 32 % | **56 %** |

The land is split sharply along the canal–desert border, which runs roughly from Chhor down past Umerkot to south-east of Kunri. Map: `04_suitability_classes.png`.

### Change since 2000

- Overall, the land got **better, not worse**:
  - The farmed zone's average score rose from **0.50 (2000–02) to 0.67 (2023–25)**.
  - About **52 % of farmed land shows a significant improvement**, probably from more irrigation and cropping.
  - About **56 % of the improving pixels still improve after removing the rainfall effect**. So part of the "greening" is just the wetter years since 2019.
- **Rainfall drives a lot of the year-to-year swings.** Farm-zone suitability follows rainfall closely (correlation 0.69). The 2002 and 2018 droughts show up as clear dips (`08_time_series.png`).
- **Fast decline is rare and patchy.** It covers 1.8 % of the land, about **110 km²** in total:
  - about 39 km² inside the farmed zone
  - about 71 km² in the desert
- When the test is re-run on 2000–2021 only, 85 % of the fast-decline pixels still show a significant decline. So these spots are not an artefact of extra satellite images being available after 2021.

### Where suitability is falling fastest (`07_decline_hotspots.png`, `decline_hotspots_5km.csv`)

1. **A north–south strip about 69.45–69.55 °E, from ~24.9 to 25.25 °N.** This is west and south-west of Kunri and south of Samaro. Five of the top eight 5 × 5 km blocks are here, with 7–11 % of each block in fast decline. Small patches inside otherwise improving farmland are typical of field-scale salinity or waterlogging, or of fields being abandoned. The data can't tell these causes apart.
2. **About 10 km east of Umerkot**, on the desert edge (block centre 25.33 °N, 69.83 °E). Only 23 % of this block is farmed. The decline here may be edge land going out of use, or the town growing, rather than salt.

### Salinity (`05_salinity_proxy.png`)

- Inside the farmed zone, the **strongest salt signal is along the eastern edge**, where irrigated land meets the desert, especially south and south-east of Kunri. Seepage and waterlogging are known problems at canal margins like this.
- The farm-zone salinity signal **jumps around a lot from year to year with no clear long-term trend**. It was highest in 2018 and 2015 and lowest in 2023.

## 4. Limitations (important for any report or presentation)

1. **No ground-truth soil samples.** "Salinity" here is a *proxy* (a satellite brightness index), not a measured EC value. To report real salinity you'd need soil EC samples, even 30–50 points. With those, the index-based score could be swapped for a trained model such as a Random Forest.
2. **Sand looks like salt.** Bright dune sand gives high SI and BI values, so the salinity map is only meaningful inside the farmed zone. The desert is greyed out.
3. **More satellites after 2021.** Landsat 9 roughly doubled the number of clear images per year from 2022. That can push "peak NDVI" up a little. The 2000–2021 re-check shows the main findings still hold.
4. **The weights are expert judgement.** Different AHP weights would shift the class boundaries. That's normal, but say so.
5. **100 m pixels** are bigger than many fields, so very small patches get blurred.
6. Town locations on the maps are approximate and only there to help you find your way.

## 5. What's in this folder

```
gee/01_extract_data.js        Earth Engine script (Step 1)
analysis/run_analysis.py      Model + trend analysis (Steps 2–3)
data/*.npy                    The 9 layers downloaded from Earth Engine
outputs/figures/              8 maps & charts (PNG)
outputs/gis/                  Results as rasters (.asc + .prj) — open in QGIS/ArcGIS
outputs/summary_by_zone.csv   Class & trend percentages per zone
outputs/decline_hotspots_5km.csv  Every 5×5 km block ranked by fast decline
outputs/yearly_means.csv      Year-by-year averages (suitability, salinity, NDVI, rain, temp)
outputs/model_info.json       AHP weights, thresholds, robustness checks
```

## 6. How to re-run it

1. **Data:** open https://code.earthengine.google.com, paste `gee/01_extract_data.js`, press **Run**, then click each `DOWNLOAD` link in the Console. Put the `.npy` files in `data/`.
2. **Analysis:** install Python 3, then run:
   ```
   pip install numpy scipy matplotlib
   python analysis/run_analysis.py data outputs
   ```
3. **Experiments to try while learning:**
   - change the AHP table (`P = ...`) to give salinity more or less weight
   - change the year range
   - move the box (`WEST/EAST/SOUTH/NORTH` in both scripts) to Badin or Mirpurkhas

## 7. Glossary

- **Remote sensing:** measuring the Earth from satellites without touching it.
- **Landsat:** NASA/USGS satellites photographing every spot every 16 days since the 1980s, at 30 m detail.
- **NDVI:** (Near-infrared − Red) / (Near-infrared + Red). Healthy plants reflect lots of near-infrared light, so higher NDVI means greener. Bare soil is about 0.1; good crops are 0.6 or more.
- **Salinity Index (SI):** √(Blue × Red). Salt crusts are bright in visible light, so a higher SI suggests more salt on bare ground.
- **Brightness Index (BI):** √(Red² + NIR²), the overall brightness of the surface.
- **Dry season (Apr–Jun):** after the wheat harvest and before the monsoon. Fields are bare and salt is easiest to see.
- **CHIRPS:** a rainfall dataset that combines satellites and rain gauges.
- **MODIS LST:** daily land surface temperature from NASA's Terra satellite.
- **AHP (Analytic Hierarchy Process):** you compare factors in pairs ("salinity is 2× as important as NDVI") and it works out the weights.
- **FAO suitability classes:** the UN Food and Agriculture Organization's standard S1/S2/S3/N land rating.
- **Mann–Kendall test:** checks whether a time series is going consistently up or down. p < 0.05 means there's under a 5 % chance it's random.
- **Sen's slope:** the median rate of change between all pairs of years. It isn't thrown off much by one odd year.
- **Pixel:** one grid square, here 0.001° ≈ 100 m.

## Data credits
- Landsat Collection 2: USGS/NASA
- CHIRPS: UCSB Climate Hazards Center
- MODIS MOD11A2: NASA LP DAAC
- SRTM: NASA/USGS
- ESA WorldCover 2021
- OpenLandMap soil texture: Hengl et al.
- JRC Global Surface Water

Salinity indices follow Khan et al. (2001, 2005) for the Indus plain. The TM/ETM+ to OLI harmonisation uses Roy et al. (2016).
