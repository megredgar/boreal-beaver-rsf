/*******************************************************************************
 * Beaver RSF (use–availability) — Alberta
 * 1) Builds a NESTED pool of available points (uniform random in area, inside
 *    the habitat mask, > EXCLUDE_RADIUS from any used point)
 * 2) Extracts distance covariates at used + available points
 * 3) Exports one CSV to Drive for model fitting in R
 *
 * The available pool is exported in random order with a 'rand' column, so any
 * N you want to test is just: avail[rank(rand) <= N]. Subsets are nested
 * (the 5k set contains the 2.5k set, etc.), which is what you want for a
 * Northrup et al. (2013)-style convergence check.
 *
 * LAND COVER: NALCMS 2020 (30 m, Landsat), already in the GEE catalog —
 *   no upload needed. ee.Image('USGS/NLCD_RELEASES/2020_REL/NALCMS')
 *   Used for: conifer, developed, and the availability mask. (Its wetland class
 *   misses treed bogs/fens, so wetlands now come from ABMI — see below.)
 *
 * BROADLEAF: SCANFI v1.2 (30 m, 2020), also in the GEE catalog.
 *   'prcB' = % broadleaf tree cover (mostly aspen/balsam poplar = beaver food).
 *   Used for: dist_deciduous (distance to pixels >= BROADLEAF_MIN % broadleaf)
 *   and pct_broadleaf (continuous % at the point).
 *   dist_deciduous and pct_broadleaf will be correlated — check |r| before
 *   putting both in one model.
 *
 * WETLAND + WATER: ABMI Alberta Wetland Inventory (10 m, Sentinel-1/2,
 *   province-wide: forested, prairie-pothole and mountain mapping areas).
 *   Codes: 0 Open water, 1 Fen, 2 Bog, 3 Marsh, 4 Swamp, 5 Upland.
 *   Used for: dist_water (ABMI open water OR NALCMS water), dist_wetland
 *   (any of fen/bog/marsh/swamp), and the wetland-type split dist_fen,
 *   dist_bog, dist_marsh, dist_swamp — all valid everywhere in Alberta.
 *   The 10 m wetland masks are aggregated to 30 m with MAX, so a small pothole
 *   anywhere in a 30 m cell still counts.
 *   dist_wetland and the four split covariates are nested — don't fit them
 *   together; the |r| check will flag this anyway.
 *
 * NICOLE !!! BEFORE RUNNING — upload these as assets (Assets tab > NEW):
 *   - beavers_AB_clean.csv (already done) -> Table asset,
 *     X = longitude, Y = latitude.
 *   - Roads: zip the .shp/.shx/.dbf/.prj together -> Table asset.
 *   - ABMI wetland raster -> Image asset. Download the .tif from
 *     https://abmi.ca/data-portal/40.html (ABMIwetlandInventory.tif, ~226 MB).
 *     On upload: Pyramiding policy = MODE, and make sure 0 is NOT set as the
 *     no-data value (0 = open water in this product!).
 ******************************************************************************/


// ============================== 0. CONFIG ==================================
var LC_ASSET    = 'USGS/NLCD_RELEASES/2020_REL/NALCMS';            // public catalog
var SCANFI_ASSET = 'projects/gcpm041u-lemur/assets/scanfi_v12/SCANFI_v1_2'; // public catalog
var WET_ASSET   = 'projects/ee-mredgar/assets/ABMIwetlandInventory';  // TODO - Nicole check path after upload
var WET_SCALE   = 10;         // ABMI wetland inventory native resolution (m)
var BROADLEAF_MIN = 25;      // % broadleaf cover for a pixel to count as "deciduous".
                              // Arbitrary — justify from lit or test 10/25/50 as sensitivity
var USED_ASSET  = 'projects/ee-mredgar/assets/beavers_AB_clean';        // TODO - Nicole change this
var ROADS_ASSET = 'projects/ee-mredgar/assets/AC_1M_Roads';       // TODO - Nicole change this 
var EXPORT_FOLDER = 'beaver_RSF';

var CRS        = 'EPSG:3400'; // Alberta 10-TM (metres, low distortion province-wide)
var LC_SCALE   = 30;          // NALCMS native resolution (m)
var DIST_SCALE = 30;          // distance surfaces computed at native resolution
var MAX_DIST   = 10000;       // distances truncated here (m). Raise if needed, but
                              // larger = slower, and >10 km is ecologically flat anyway
var N_AVAIL    = 50000;       // size of the full available pool (test subsets in R)
var OVERSAMPLE = 1.6;         // extra candidates to survive mask + exclusion filtering
var EXCLUDE_RADIUS = 100;     // m; no available points within this of a used point
var SEED = 2026;

// ---- NALCMS 2020 class codes (confirm against the catalog class table) ----
//  1 Temperate/sub-polar needleleaf forest   2 Sub-polar taiga needleleaf forest
//  5 Temperate/sub-polar broadleaf deciduous 6 Mixed forest
//  8 Temperate/sub-polar shrubland          10 Temperate/sub-polar grassland
// 14 Wetland  15 Cropland  16 Barren  17 Urban/built-up  18 Water  19 Snow/ice
var LC = {
  water:      [18],
  developed:  [17],
  deciduous:  [5],
  coniferous: [1, 2],
  mixedwood:  [6],
  wetland:    [14],
  shrub:      [8],
  agriculture:[15]
};
// ---- ABMI Alberta Wetland Inventory class codes (from the ABMI metadata) ----
// 0 Open water  1 Fen  2 Bog  3 Marsh  4 Swamp  5 Upland
var WET = {
  water:  [0],
  fen:    [1],
  bog:    [2],
  marsh:  [3],
  swamp:  [4],
  anyWet: [1, 2, 3, 4]
};

// Classes that count as "available" habitat (per the proposal's land mask)
var AVAILABLE_CLASSES = [].concat(
  LC.water, LC.wetland, LC.shrub, LC.coniferous,
  LC.deciduous, LC.mixedwood, LC.agriculture
);


// ============================ 1. INPUTS ====================================
var alberta = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM1_NAME', 'Alberta'));       // or upload your own boundary
var region = alberta.geometry();

var lc    = ee.Image(LC_ASSET).select('landcover').rename('lc');
var prcB  = ee.Image(SCANFI_ASSET).select('prcB').unmask(0);   // % broadleaf cover
var wet   = ee.Image(WET_ASSET).select(0).rename('wet');       // ABMI wetland inventory

var roads = ee.FeatureCollection(ROADS_ASSET);
var used  = ee.FeatureCollection(USED_ASSET)
  .filterBounds(region)
  .map(function (f) { return f.set('used', 1, 'rand', -1); });

print('Used points (after Alberta filter):', used.size());

// Sanity check — which NALCMS codes occur around Edmonton (should include
// 1/5/6/14/15/17/18). Run once, then comment out if you like.
var testBox = ee.Geometry.Rectangle([-114.5, 53.0, -112.5, 54.5]);
print('LC class histogram (test box):',
  lc.reduceRegion({reducer: ee.Reducer.frequencyHistogram(),
                   geometry: testBox, scale: 100, maxPixels: 1e9}));

// Same check for the ABMI wetland raster — should show keys 0–5.
// If '0' (open water) is missing, the upload treated 0 as no-data: re-upload.
print('ABMI wetland class histogram (test box):',
  wet.reduceRegion({reducer: ee.Reducer.frequencyHistogram(),
                    geometry: testBox, scale: 100, maxPixels: 1e9}));


// ========================= 2. HELPER FUNCTIONS =============================
// 0/1 image for a set of class codes
function classMask(codes) {
  return lc.remap(codes, ee.List.repeat(1, codes.length), 0);
}

// Euclidean distance (m) to the nearest 1-pixel of any 0/1 image, at 30 m.
// Truncated at MAX_DIST.
function distToMask(mask, name) {
  var proj = ee.Projection(CRS).atScale(DIST_SCALE);
  var nPix = Math.ceil(MAX_DIST / DIST_SCALE);
  return mask.reproject(proj)
    .fastDistanceTransform(nPix, 'pixels', 'squared_euclidean')
    .sqrt().multiply(DIST_SCALE)
    .reproject(proj)
    .unmask(MAX_DIST).min(MAX_DIST)
    .rename(name);
}

// Distance to a set of NALCMS classes
function distToClass(codes, name) {
  return distToMask(classMask(codes), name);
}

// 0/1 image for a set of ABMI wetland classes, aggregated 10 m -> 30 m with
// MAX (a pothole anywhere in the 30 m cell counts).
function wetMask(codes) {
  return wet.remap(codes, ee.List.repeat(1, codes.length), 0)
    .unmask(0)
    .reproject(ee.Projection(CRS).atScale(WET_SCALE))
    .reduceResolution({reducer: ee.Reducer.max(), maxPixels: 16});
}


// ======================== 3. COVARIATE STACK ===============================
var distRoads = roads.distance({searchRadius: MAX_DIST, maxError: 10})
  .unmask(MAX_DIST).min(MAX_DIST).rename('dist_roads');

var covs = ee.Image.cat([
  distToMask(wetMask(WET.water).or(classMask(LC.water)), 'dist_water'),  // ABMI + NALCMS
  distToMask(wetMask(WET.anyWet), 'dist_wetland'),   // ABMI fen/bog/marsh/swamp
  distToMask(wetMask(WET.fen),    'dist_fen'),
  distToMask(wetMask(WET.bog),    'dist_bog'),
  distToMask(wetMask(WET.marsh),  'dist_marsh'),
  distToMask(wetMask(WET.swamp),  'dist_swamp'),
  wet.rename('wet_class'),                  // ABMI class at the point (QA)
  distToMask(prcB.gte(BROADLEAF_MIN), 'dist_deciduous'),   // SCANFI broadleaf
  distToClass(LC.coniferous, 'dist_coniferous'),
  distToClass(LC.developed,  'dist_developed'),
  distRoads,
  prcB.rename('pct_broadleaf'),             // SCANFI % broadleaf at the point
  lc.rename('lc_class')                     // handy for QA / summaries
]);
var COV_NAMES = ['dist_water', 'dist_wetland', 'dist_fen', 'dist_bog',
                 'dist_marsh', 'dist_swamp', 'wet_class', 'dist_deciduous',
                 'dist_coniferous', 'dist_developed', 'dist_roads',
                 'pct_broadleaf', 'lc_class'];


// ===================== 4. AVAILABLE POINT POOL =============================
// Habitat mask, minus a buffer around every used point
var nearUsed = ee.Image(0).byte()
  .paint(used.map(function (f) { return f.buffer(EXCLUDE_RADIUS); }), 1);

var availMask = classMask(AVAILABLE_CLASSES)
  .and(nearUsed.not())
  .selfMask()
  .rename('keep');

// Uniform-in-area candidates -> keep those that land in the mask
var candidates = ee.FeatureCollection.randomPoints({
  region: region,
  points: Math.round(N_AVAIL * OVERSAMPLE),
  seed: SEED,
  maxError: 10
});

var avail = availMask.sampleRegions({
    collection: candidates, scale: LC_SCALE, projection: CRS,
    geometries: true, tileScale: 16
  })                                  // masked (non-habitat / near-used) points drop out
  .randomColumn('rand', SEED)
  .sort('rand')
  .limit(N_AVAIL)
  .map(function (f) { return f.set('used', 0, 'id', -1); });

// Should be 50,000. If it's short, raise OVERSAMPLE (e.g. to 2) and rerun.
print('Available n:', avail.size());


// ====================== 5. EXTRACT + EXPORT ================================
var allPts = used.merge(avail).map(function (f) {
  var xy = f.geometry().coordinates();
  return f.set('lon', xy.get(0), 'lat', xy.get(1));
});

var out = covs.sampleRegions({
  collection: allPts,
  properties: ['id', 'used', 'rand', 'lon', 'lat'],
  scale: LC_SCALE,
  projection: CRS,
  geometries: false,
  tileScale: 16
});

Export.table.toDrive({
  collection: out,
  description: 'beaver_RSF_used_avail_covariates',
  folder: EXPORT_FOLDER,
  fileFormat: 'CSV',
  selectors: ['id', 'used', 'rand', 'lon', 'lat'].concat(COV_NAMES)
});

// Optional: export the available points alone (check them on a map in QGIS)
Export.table.toDrive({
  collection: avail,
  description: 'beaver_RSF_available_points',
  folder: EXPORT_FOLDER,
  fileFormat: 'SHP'
});


// ============================ 6. QA MAP ====================================
// Don't print() the big collections — it will time out. Look, then run Tasks.
Map.centerObject(testBox, 9);
Map.addLayer(lc, {min: 1, max: 19, palette: [
  '033e00','939b71','196d12','1fab01','5b725c','6b7d2c','b29d29','b48833',
  'e9da5d','e0cd88','a07451','bad292','3f8970','6ca289','e6ad6a','a9abae',
  'db2126','4c73a1','fff7fe']}, 'NALCMS 2020', false);
Map.addLayer(covs.select('dist_water'), {min: 0, max: 2000}, 'dist_water');
Map.addLayer(wet, {min: 0, max: 5, palette: [
  '0070FF','E6A800','4C0073','00A884','5C8944','FFFFFF']},   // water fen bog marsh swamp upland
  'ABMI wetland inventory', false);
Map.addLayer(covs.select('dist_wetland'), {min: 0, max: 2000}, 'dist_wetland', false);
Map.addLayer(prcB.selfMask(), {min: 0, max: 100, palette: ['yellow', 'darkgreen']},
             'SCANFI % broadleaf', false);
Map.addLayer(covs.select('dist_deciduous'), {min: 0, max: 2000}, 'dist_deciduous', false);
Map.addLayer(distRoads, {min: 0, max: 5000}, 'dist_roads', false);
Map.addLayer(used, {color: 'red'}, 'used');
