import geopandas as gpd
import shapely
import requests
import rasterio
import os

# check if ./temp exists
if not os.path.exists('temp'):
    os.mkdir('temp')

file_name = 'mydata/border_2180.shp'

df = gpd.read_file(file_name)
# reproject to epsg:2180
df = df.to_crs(epsg=2180)

print(df)

# make a buffer of feature with name Świeradów Zdrój
feature = df[df['JPT_NAZWA_'] == 'Świeradów-Zdrój']
buffer = feature.buffer(200)

# save to /temp/buffer.shp
buffer.to_file('temp/buffer.shp')

buffer_bbox = buffer.total_bounds
center = buffer.centroid


# download DEM data from wcs

version = '1.0.0'
request_type = 'GetCoverage'
format = 'image/tiff'
coverage = 'DTM_PL-KRON86-NH_TIFF'
bbox_str = str(buffer_bbox).replace(' ', ',')[1:-1]
crs = 'EPSG:2180'
response_crs = 'EPSG:2180'

params = {
    'SERVICE': 'WCS',
    'VERSION': version,
    'REQUEST': request_type,
    'FORMAT': format,
    'COVERAGE': coverage,
    'BBOX': bbox_str,
    'CRS': crs,
    'RESPONSE_CRS': response_crs,
    'WIDTH': (buffer_bbox[2] - buffer_bbox[0])//5,
    'HEIGHT': (buffer_bbox[3] - buffer_bbox[1])//5
}

url = 'https://mapy.geoportal.gov.pl/wss/service/PZGIK/NMT/GRID1/WCS/DigitalTerrainModelFormatTIFF?'

for k, v in params.items():
    url += f'{k}={v}&'
url = url[:-1]

response = requests.get(url)
data = response.content


file_path = 'temp/xd.tif'

with open(file_path, 'wb') as f:
    f.write(data)

luban_gml = '0210_GML'
obok_gml = '0212_GML'


