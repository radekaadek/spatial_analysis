from typing import Sequence
import geopandas as gpd
import rasterio
import os
import pandas as pd
import numpy as np

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

def process_feature(feature_name: str, area_polygon: gpd.GeoDataFrame) -> pd.DataFrame:
    luban_name = '02_GML/0210_GML'
    lwowecki_name = '02_GML/0212_GML'
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


area_polygon = gpd.read_file('swieradow_buffer/swieradow_buffer.shp')

rivers_frame = process_feature(rivers, area_polygon)
buildings_frame = process_feature(buildings, area_polygon)
woda_powierzchnia_frame = process_feature(woda_powierzchnia, area_polygon)
zabudowa_frame = process_feature(zabudowa, area_polygon)
teren_lesny_frame = process_feature(teren_lesny, area_polygon)
roslinnosc_krzewiasta_frame = process_feature(roslinnosc_krzewiasta, area_polygon)
teren_pod_drogami_frame = process_feature(teren_pod_drogami, area_polygon)
grunt_nieuzytkowy_frame = process_feature(grunt_nieuzytkowy, area_polygon)
plac_frame = process_feature(plac, area_polygon)
skladowisko_odpadow_frame = process_feature(skladowisko_odpadow, area_polygon)
wyrobisko_frame = process_feature(wyrobisko, area_polygon)


# load the dem
result = rasterio.open('temp/xd.tif')
band = result.read(1)
band = np.zeros(band.shape)

# save buildings to buildings.gpkg
# buildings_frame.to_file('buildings.gpkg')

