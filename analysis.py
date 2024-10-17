from typing import Sequence
import geopandas as gpd
import rasterio
import os
import pandas as pd

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
    # return feature.clip(area_polygon)
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

area_polygon = gpd.read_file('swieradow_buffer/swieradow_buffer.shp')

rivers_frame = process_feature(rivers, area_polygon)
buildings_frame = process_feature(buildings, area_polygon)

print(rivers_frame.head())
print(buildings_frame.head())

# save buildings to buildings.gpkg
buildings_frame.to_file('buildings.gpkg')

