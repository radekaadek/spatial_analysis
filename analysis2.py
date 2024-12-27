import math
import subprocess
import os
import shapely
import rasterio
import rasterstats
import numpy as np
import pandas as pd
import geopandas as gpd
from typing import Sequence

def find_file(base_path: str, file_ending: str) -> str:
    for root, _, files in os.walk(base_path):
        for file in files:
            if file.endswith(file_ending):
                return os.path.join(root, file)
    raise FileNotFoundError(f'File with ending {file_ending} not found')

def read_and_concat_files(file_names: Sequence[str]) -> pd.DataFrame:
    gdfs = [gpd.read_file(file_name) for file_name in file_names]
    return pd.concat(gdfs)

# get area of analysis

luban_name = '02_GML/0210_GML'
lwowecki_name = '02_GML/0212_GML'
def process_feature(feature_name: str, area_polygon: gpd.GeoDataFrame) -> pd.DataFrame:
    f1 = find_file(luban_name, feature_name)
    f2 = find_file(lwowecki_name, feature_name)
    feature = read_and_concat_files([f1, f2])
    # crop feature to area of interest
    clipped = feature.clip(area_polygon)
    if clipped is None:
        return gpd.GeoDataFrame()
    return clipped





# base_name = 'PL.PZGiK.337.'
# luban = '0210_GML'
# lwowecki = '0212_GML'

# get rivers
luban_name = '02_GML/0210_GML'
lwowecki_name = '02_GML/0212_GML'

rivers = 'SWRS_L.xml'
buildings = 'BUBD_A.xml'

woda_powierzchnia = 'PTWP_A.xml'
zabudowa = 'PTZB_A.xml'
teren_lesny = 'PTLZ_A.xml'
roslinnosc_krzewiasta = 'PTRK_A.xml'
uprawa_trwala = 'PTUT_A.xml'
roslinnosc_trawiasta = 'PTTR_A.xml'
teren_pod_drogami = 'PTKM_A.xml'
grunt_nieuzytkowy = 'PTGN_A.xml'
plac = 'PTPL_A.xml'
skladowisko_odpadow = 'PTSO_A.xml'
wyrobisko = 'PTWZ_A.xml'
pozostale = 'PTNZ_A.xml'
droga = 'SKDR_L.xml'
power_line = 'SULN_L.xml'


name_to_xml_name = {
    'rivers': 'SWRS_L.xml',
    'buildings': 'BUBD_A.xml',
    'woda_powierzchnia': 'PTWP_A.xml',
    'zabudowa': 'PTZB_A.xml',
    'teren_lesny': 'PTLZ_A.xml',
    'roslinnosc_krzewiasta': 'PTRK_A.xml',
    'teren_pod_drogami': 'PTKM_A.xml',
    'grunt_nieuzytkowy': 'PTGN_A.xml',
    'plac': 'PTPL_A.xml',
    'skladowisko_odpadow': 'PTSO_A.xml',
    'wyrobisko': 'PTWZ_A.xml',
    'pozostale': 'PTNZ_A.xml',
    'power_line': 'SULN_L.xml'
}

vector_dir = 'vectors'

# area_polygon = gpd.read_file('swieradow_buffer/swieradow_buffer.shp')
# load the dem
dem_path = 'rasters/dem_original.tif'
result = rasterio.open(dem_path)
band = result.read(1)
band_shape = band.shape
bounds = result.bounds
area_polygon = gpd.GeoDataFrame(geometry=[shapely.geometry.box(*bounds)]).set_crs("EPSG:2180")

rivers_frame = process_feature(rivers, area_polygon)
buildings_frame = process_feature(buildings, area_polygon)
# select funkcjaOgolnaBudynku == 'budynki mieszkalne'
buildings_frame = buildings_frame[buildings_frame['funkcjaOgolnaBudynku'] == 'budynki mieszkalne']
woda_powierzchnia_frame = process_feature(woda_powierzchnia, area_polygon)
zabudowa_frame = process_feature(zabudowa, area_polygon)
teren_lesny_frame = process_feature(teren_lesny, area_polygon)
roslinnosc_krzewiasta_frame = process_feature(roslinnosc_krzewiasta, area_polygon)
teren_pod_drogami_frame = process_feature(teren_pod_drogami, area_polygon)
grunt_nieuzytkowy_frame = process_feature(grunt_nieuzytkowy, area_polygon)
plac_frame = process_feature(plac, area_polygon)
skladowisko_odpadow_frame = process_feature(skladowisko_odpadow, area_polygon)
wyrobisko_frame = process_feature(wyrobisko, area_polygon)
droga_frame = process_feature(droga, area_polygon)
power_line_frame = process_feature(power_line, area_polygon)
# read vectors/dzialki.gpkg
dzialki = gpd.read_file(f"{vector_dir}/dzialki.gpkg")
dzialki_frame = dzialki.clip(area_polygon)


all_roads = droga_frame.copy()

good_roads = droga_frame
good_roads = good_roads[good_roads['kategoriaZarzadzania'].isin(['wojewódzka'])]
roads_geometry = droga_frame['geometry']
roads_now = gpd.GeoDataFrame(columns=['geometry'], geometry='geometry', crs=droga_frame.crs)
intersections = gpd.GeoDataFrame(columns=['geometry'], geometry='geometry', crs=droga_frame.crs)
for i, row in good_roads.iterrows():
    inter = row['geometry'].intersection(roads_geometry)
    # delete all geometries that are not points
    inter = inter[~inter.is_empty]
    inter = inter[inter.geom_type == 'Point']
    try: 
        g = next(iter(inter))
        intersections.at[i, 'geometry'] = g
    except StopIteration:
        pass
    roads_now._append(row)


drogi_hard_names = ['masa bitumiczna', 'beton', 'kostka prefabrykowana', 'tłuczeń', 'płyty betonowe', 'kostka kamienna', 'bruk']
for name in drogi_hard_names:
    droga_frame = droga_frame[droga_frame['materialNawierzchni'] != name]

# saave all frames to file
if not os.path.exists(vector_dir):
    os.mkdir(vector_dir)

rivers_frame.to_file(f'{vector_dir}/rivers.gpkg')
buildings_frame.to_file(f'{vector_dir}/buildings.gpkg')
woda_powierzchnia_frame.to_file(f'{vector_dir}/woda_powierzchnia.gpkg')
zabudowa_frame.to_file(f'{vector_dir}/zabudowa.gpkg')
teren_lesny_frame.to_file(f'{vector_dir}/teren_lesny.gpkg')
roslinnosc_krzewiasta_frame.to_file(f'{vector_dir}/roslinnosc_krzewiasta.gpkg')
teren_pod_drogami_frame.to_file(f'{vector_dir}/teren_pod_drogami.gpkg')
grunt_nieuzytkowy_frame.to_file(f'{vector_dir}/grunt_nieuzytkowy.gpkg')
plac_frame.to_file(f'{vector_dir}/plac.gpkg')
skladowisko_odpadow_frame.to_file(f'{vector_dir}/skladowisko_odpadow.gpkg')
wyrobisko_frame.to_file(f'{vector_dir}/wyrobisko.gpkg')
droga_frame.to_file(f'{vector_dir}/droga.gpkg')
all_roads.to_file(f'{vector_dir}/all_roads.gpkg')
intersections.to_file(f'{vector_dir}/intersections.gpkg')
dzialki_frame.to_file(f'{vector_dir}/dzialki.gpkg')
power_line_frame.to_file(f'{vector_dir}/power_line.gpkg')

#### PARAMETERS ####
pixel_size = 5
weights: dict[str, float] = {
    'dzialki': 1,
    'odleglosc_budynkow': 1,
    'woda': 1,
    'las': 6,
    'nachylenie': 1,
    'wystawa': 1,
    'dostep_drog': 1,
    'skrzyzowania': 1,
}
# select the top % of all values, including nan values
percent_value = 0.65
# procent of pixles inside parcels
procent_of_pixels = 0.60
dzialka_area_min = 20000



# normalize weights so that sum of weights is 1
weights_sum = sum(weights.values())
for key, value in weights.items():
    weights[key] = value/weights_sum


raster_dir = 'rasters'
if not os.path.exists(raster_dir):
    os.mkdir(raster_dir)

# rasterize all features
for vector_file in os.listdir(vector_dir):
    file_path = f'{vector_dir}/{vector_file}'
    feature = vector_file.split('.')[0]
    gdal_command = f"gdal_rasterize -burn 1 -ts {band_shape[1]} {band_shape[0]} -te {bounds[0]} {bounds[1]} {bounds[2]} {bounds[3]} {file_path} {raster_dir}/{feature}.tif"
    subprocess.run(gdal_command, shell=True)

def calculate_distance(raster_path: str, output_path: str) -> None:
    gdal_command = f"gdal_proximity.py {raster_path} {output_path}"
    subprocess.run(gdal_command, shell=True)

distances_dir = 'distances'

if not os.path.exists(distances_dir):
    os.mkdir(distances_dir)

for file in os.listdir(raster_dir):
    file_path = f'{raster_dir}/{file}'
    feature = file.split('.')[0]
    output_path = f'{distances_dir}/{feature}.tif'
    calculate_distance(file_path, output_path)

# dem slope raster
dem_raster_path = f'{raster_dir}/dem_original.tif'
gdal_command = f'gdaldem slope {dem_raster_path} {raster_dir}/slope.tif'
os.system(gdal_command)

# dem aspect raster
gdal_command = f'gdaldem aspect {dem_raster_path} {raster_dir}/aspect.tif'
os.system(gdal_command)

result_band = np.zeros(band.shape)

#### BUILDINGS ####
# buildings - at <150 meters set nodata
buildings_raster = rasterio.open(f'{distances_dir}/buildings.tif')
building_distances = buildings_raster.read(1)*pixel_size
result_band[building_distances <= 150] = np.nan
penalty_mask = (150 < building_distances) & (building_distances < 200)
result_band[penalty_mask] += (400 - building_distances[penalty_mask] * 2)*weights['dzialki']

#### WATER ####
# water - at <100 meters set nodata
# the closer to the water, the better
rivers_raster = rasterio.open(f'{distances_dir}/rivers.tif')
water_distances = rivers_raster.read(1)*pixel_size
lake_raster = rasterio.open(f'{distances_dir}/woda_powierzchnia.tif')
water_distances = np.minimum(water_distances, lake_raster.read(1)*pixel_size)
result_band[water_distances <= 100] = np.nan
# dont forget to check if nodata is set
penalty_mask = (~np.isnan(result_band))
result_band[penalty_mask] += (water_distances[penalty_mask] * 0.5)*weights['woda']

#### BORDER ####
border_raster = rasterio.open(f'{distances_dir}/border_2180.tif')
result_band[border_raster.read(1) > 0] = np.nan

#### FOREST ####
forest_raster = rasterio.open(f'{distances_dir}/teren_lesny.tif')
# forest at <15 meters set nodata
forest_distances = forest_raster.read(1)*pixel_size
result_band[forest_distances <= 15] = np.nan
# forest at >15 <100
max_forest = 100
penalty_mask = (forest_distances > 15) & (forest_distances < max_forest)
result_band[penalty_mask] += (max_forest - forest_distances[penalty_mask])*weights['las']


#### ROADS ####
roads_raster = rasterio.open(f'{distances_dir}/all_roads.tif')
roads_distances = roads_raster.read(1)*pixel_size
road_mask = (roads_distances <= 10)
result_band[road_mask] = np.nan
penalty_mask = (roads_distances > 15)
result_band[penalty_mask] += (roads_distances[penalty_mask] * 0.5)*weights['dostep_drog']


#### SLOPE ####
slope_raster = rasterio.open(f'{distances_dir}/slope.tif')
slope_deg = slope_raster.read(1)
result_band[slope_deg > 10] = np.nan
penalty_mask = (slope_deg <= 10)
result_band[penalty_mask] += (slope_deg[penalty_mask] * 4)*weights['nachylenie']

#### ASPECT ####
aspect_raster = rasterio.open(f'{raster_dir}/aspect.tif')
aspect_deg = aspect_raster.read(1)
# optymalnie: stoki poludniowe SW-SE
result_band[aspect_deg > 0] += (abs(aspect_deg[aspect_deg > 0] - 180))*weights['wystawa']

#### INTERSECTIONS ####
# the further away from an intersection the worse
intersections_raster = rasterio.open(f'{distances_dir}/intersections.tif')
intersections_distances = intersections_raster.read(1)*pixel_size
penalty_mask = (intersections_distances >= 0)
result_band[penalty_mask] += (intersections_distances[penalty_mask] * 0.25)*weights['skrzyzowania']


# minmax scale
minmax_scaled = np.nanmin(result_band), np.nanmax(result_band)
result_band = (result_band - minmax_scaled[0]) / (minmax_scaled[1] - minmax_scaled[0])
result_band = 1 - result_band

result_band = np.nan_to_num(result_band)

result_profile = result.profile
with rasterio.open("result.tif", "w", **result_profile) as dst:
    dst.write(result_band, 1)


# set to 1 or 0
result_band[result_band >= percent_value] = 1
result_band[result_band < percent_value] = 0

with rasterio.open("result2.tif", "w", **result_profile) as dst:
    dst.write(result_band, 1)

# sieve pixels
min_area = 50
pixel_area = pixel_size * pixel_size
pixels = math.ceil(min_area / pixel_area)
gdal_command = f"gdal_sieve.py -st 5 result2.tif result3.tif"
subprocess.run(gdal_command, shell=True)

# read result3 with rasterio
result3 = rasterio.open("result3.tif")
result3_band = result3.read(1)
t = result3.transform

# get indices where result3 is 1
result3_indexes = np.where(result_band == 1) # results before sieving

points = list(zip(*rasterio.transform.xy(t, result3_indexes[0], result3_indexes[1])))

xs, ys = zip(*points)
shapely_points = [shapely.geometry.Point(x, y) for x, y in zip(xs, ys)]
points_series = gpd.GeoSeries(shapely_points, crs=result3.crs)

points_gdf = gpd.GeoDataFrame(geometry=points_series)

# joint points with buildings
dzialki_with_points = dzialki_frame.sjoin(points_gdf, how='inner', predicate='intersects')

# calcualte raster stats for dzialki_with_points.gpkg
# get affine transform from result2.tif
with rasterio.open("result2.tif") as src:
    affine = src.transform
stats = rasterstats.zonal_stats(dzialki_with_points, result_band, affine=affine, stats=['mean'])

# add to dzialki_with_points.gpkg
for (i, row), stats in zip(dzialki_with_points.iterrows(), stats):
    dzialki_with_points.at[i, 'mean'] = stats['mean']

# chooose dzialki with mean > 0.60
dzialki_with_points = dzialki_with_points[dzialki_with_points['mean'] > procent_of_pixels]
# take działki with area > 50
dzialki_with_points['area'] = dzialki_with_points.geometry.area
# delete dzialki with duplicate ID_DZIALKI
dzialki_with_points = dzialki_with_points.drop_duplicates(subset='ID_DZIALKI')
dzialki_with_points.to_file('dzialki_mean.gpkg')


# load kody.geojson
kody = gpd.read_file('kody.geojson')
unique_kody = list(kody['X_KOD'].unique())
kod_name = {i+1: kod for i, kod in enumerate(unique_kody)}
# reverse dictionary
name_kod = {v: k for k, v in kod_name.items()}

def get_kod_name(kod: str) -> int:
    x_kod_to_koszt = {
        "PTWP01": 99999,
        "PTWP02": 200,
        "PTWP03": 99999,
        "PTZB02": 100,
        "PTZB01": 200,
        "PTZB05": 50,
        "PTZB04": 200,
        "PTZB03": 200,
        "PTLZ01": 100,
        "PTLZ02": 50,
        "PTLZ03": 50,
        "PTRK01": 15,
        "PTRK02": 15,
        "PTUT03": 100,
        "PTUT02": 90,
        "PTUT04": 20,
        "PTUT05": 20,
        "PTUT01": 99999,
        "PTTR02": 1,
        "PTTR01": 20,
        "PTKM02": 200,
        "PTKM01": 100,
        "PTKM03": 200,
        "PTKM04": 99999,
        "PTGN01": 1,
        "PTGN02": 1,
        "PTGN03": 1,
        "PTGN04": 1,
        "PTPL01": 50,
        "PTSO01": 99999,
        "PTSO02": 99999,
        "PTWZ01": 99999,
        "PTWZ02": 99999,
        "PTNZ01": 150,
        "PTNZ02": 150,
    }
    return x_kod_to_koszt[kod]

kody['kod_name'] = kody['X_KOD'].apply(get_kod_name)
# set kod_name to int16
kody['kod_name'] = kody['kod_name'].astype(np.int16)

kody.to_file('kody2.geojson')

# rasterize
command = f"rio rasterize --profile dtype=int16 --property kod_name --overwrite --like result.tif kody2.geojson stuff.tif"
subprocess.run(command, shell=True)


