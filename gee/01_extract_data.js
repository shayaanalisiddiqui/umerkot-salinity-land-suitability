/**************************************************************************
 * Umerkot / Nara Canal (Sindh, Pakistan)
 * Agricultural land suitability & salinity trend project
 * STEP 1: Extract yearly satellite layers (2000–2025) with Google Earth Engine
 *
 * Paste into https://code.earthengine.google.com and press Run.
 * The Console prints one download link per layer (NumPy .npy files).
 *
 * Grid used for every layer (so all layers line up pixel-for-pixel):
 *   EPSG:4326, top-left = (69.3 E, 25.6 N), pixel = 0.001 deg (~100 m)
 *   800 columns x 700 rows  ->  69.3–70.1 E, 24.9–25.6 N
 **************************************************************************/

var WEST = 69.3, EAST = 70.1, SOUTH = 24.9, NORTH = 25.6;
var aoi = ee.Geometry.Rectangle([WEST, SOUTH, EAST, NORTH], null, false);
var START_YEAR = 2000, END_YEAR = 2025;

var GRID = {crs: 'EPSG:4326', crs_transform: [0.001, 0, WEST, 0, -0.001, NORTH], dimensions: '800x700'};
var GRID_FINE = {crs: 'EPSG:4326', crs_transform: [0.0005, 0, WEST, 0, -0.0005, NORTH], dimensions: '1600x1400'};
var NODATA = -32768;

Map.centerObject(aoi, 9);
Map.addLayer(aoi, {color: 'red'}, 'Study area', true, 0.3);

/* ---------------- Landsat 5/7/8/9 Collection 2 surface reflectance ---------------- */
// Cloud / shadow mask + scale factors
function prepL2(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)   // dilated cloud
      .and(qa.bitwiseAnd(1 << 3).eq(0))     // cloud
      .and(qa.bitwiseAnd(1 << 4).eq(0));    // cloud shadow
  var sat = img.select('QA_RADSAT').eq(0);
  return img.updateMask(clear).updateMask(sat);
}
// TM / ETM+ -> common band names, then harmonise to OLI (Roy et al. 2016)
function prepTM(img) {
  img = prepL2(img);
  var sr = img.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'])
      .multiply(0.0000275).add(-0.2);
  var slopes = ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]);
  var itcp = ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]);
  sr = sr.multiply(slopes).add(itcp)
      .rename(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']);
  return sr.copyProperties(img, ['system:time_start']);
}
function prepOLI(img) {
  img = prepL2(img);
  var sr = img.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'])
      .multiply(0.0000275).add(-0.2)
      .rename(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']);
  return sr.copyProperties(img, ['system:time_start']);
}

var L5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(aoi).map(prepTM);
// Landsat 7 only up to 2012 (before Landsat 8), to avoid SLC-off stripes & orbit drift
var L7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(aoi)
    .filterDate('2000-01-01', '2013-01-01').map(prepTM);
var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(aoi).map(prepOLI);
var L9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(aoi).map(prepOLI);
var landsat = L5.merge(L7).merge(L8).merge(L9);

function addIndices(img) {
  var ndvi = img.normalizedDifference(['nir', 'red']).rename('ndvi');
  // Salinity Index  SI = sqrt(Blue * Red)            (Khan et al. 2001/2005, Indus plain)
  var si = img.expression('sqrt(max(b,0) * max(r,0))', {b: img.select('blue'), r: img.select('red')}).rename('si');
  // Brightness Index BI = sqrt(Red^2 + NIR^2)        (salt crusts are bright)
  var bi = img.expression('sqrt(r*r + n*n)', {r: img.select('red'), n: img.select('nir')}).rename('bi');
  return img.addBands([ndvi, si, bi]).toFloat();
}
landsat = landsat.map(addIndices);

/* ---------------- Yearly layers ---------------- */
// Guarantees the band exists even if a period has zero clear scenes (gives a masked pixel)
var EMPTY = ee.Image.constant([0, 0, 0, 0, 0, 0, 0, 0, 0]).toFloat()
    .rename(['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'ndvi', 'si', 'bi']).updateMask(0);
function period(start, end) {
  return landsat.filterDate(start, end).merge(ee.ImageCollection([EMPTY]));
}

function yearStack(fn, prefix) {
  var bands = [];
  for (var y = START_YEAR; y <= END_YEAR; y++) {
    bands.push(fn(y).rename(prefix + '_' + y));
  }
  return ee.Image.cat(bands);
}

// (a) Peak greenness of the year: 90th-percentile NDVI (crop vigour / cropping)
var ndviPeak = yearStack(function (y) {
  return period(y + '-01-01', (y + 1) + '-01-01').select('ndvi')
      .reduce(ee.Reducer.percentile([90]));
}, 'ndvi_p90');

// (b) Dry-season (Apr 1 – Jun 15) medians: fields mostly bare, salt crusts most visible
function drySeason(y) {
  return period(y + '-04-01', y + '-06-16');
}
var siDry = yearStack(function (y) { return drySeason(y).select('si').median(); }, 'si_dry');
var biDry = yearStack(function (y) { return drySeason(y).select('bi').median(); }, 'bi_dry');
var ndviDry = yearStack(function (y) { return drySeason(y).select('ndvi').median(); }, 'ndvi_dry');

// (c) Annual rainfall (CHIRPS, mm)
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD');
var rain = yearStack(function (y) {
  return chirps.filterDate(y + '-01-01', (y + 1) + '-01-01').select('precipitation').sum();
}, 'rain');

// (d) Mean day-time land surface temperature, Apr–Jun (MODIS, deg C)
var modis = ee.ImageCollection('MODIS/061/MOD11A2').select('LST_Day_1km');
var lst = yearStack(function (y) {
  return modis.filterDate(y + '-04-01', y + '-07-01').mean().multiply(0.02).subtract(273.15);
}, 'lst');

// (e) Number of clear Landsat scenes per year (data-quality check)
var nScenes = yearStack(function (y) {
  return period(y + '-01-01', (y + 1) + '-01-01').select('ndvi').count().unmask(0);
}, 'n');

/* ---------------- Static layers ---------------- */
var dem = ee.Image('USGS/SRTMGL1_003');
var slope = ee.Terrain.slope(dem).rename('slope_deg');
var sand = ee.Image('OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02').select('b0').rename('sand_pct');
var clay = ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02').select('b0').rename('clay_pct');
var wc = ee.ImageCollection('ESA/WorldCover/v200').first().rename('worldcover');
var water = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0).rename('water_occ');
var staticImg = ee.Image.cat([dem.rename('elev_m'), slope, sand, clay, wc, water]).toFloat();

// True-colour image of the dry season 2024–25 for figures (~50 m)
var rgb = landsat.filter(ee.Filter.calendarRange(4, 6, 'month'))
    .filterDate('2024-01-01', '2026-01-01').median()
    .select(['red', 'green', 'blue'])
    .unitScale(0, 0.35).multiply(255).clamp(0, 255).toUint8();

/* ---------------- Preview on the map ---------------- */
Map.addLayer(rgb, {}, 'True colour 2024-25 (dry season)');
Map.addLayer(ndviPeak.select('ndvi_p90_2025'), {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green']}, 'Peak NDVI 2025', false);
Map.addLayer(siDry.select('si_dry_2025'), {min: 0.05, max: 0.3, palette: ['blue', 'white', 'red']}, 'Salinity index 2025', false);

/* ---------------- Download links ---------------- */
function link(img, name, grid) {
  var params = {name: name, format: 'NPY', crs: grid.crs,
                crs_transform: grid.crs_transform, dimensions: grid.dimensions};
  img.getDownloadURL(params, function (url, err) {
    if (err) { print('ERROR ' + name + ': ' + err); return; }
    print(ui.Label('DOWNLOAD ' + name, {color: 'blue'}, url));
  });
}
function asInt16(img, scale) {
  return img.multiply(scale).round().unmask(NODATA).toInt16();
}

link(asInt16(ndviPeak, 10000), 'ndvi_p90', GRID);
link(asInt16(ndviDry, 10000), 'ndvi_dry', GRID);
link(asInt16(siDry, 10000), 'si_dry', GRID);
link(asInt16(biDry, 10000), 'bi_dry', GRID);
link(asInt16(rain, 1), 'rain_mm', GRID);
link(asInt16(lst, 100), 'lst_c_x100', GRID);
link(asInt16(nScenes, 1), 'n_scenes', GRID);
link(staticImg.unmask(NODATA), 'static', GRID);
link(rgb.unmask(0), 'rgb_2025', GRID_FINE);
