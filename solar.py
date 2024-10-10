import geopandas as gpd
import shapely

file_name = 'mydata/border_2180.shp'

df = gpd.read_file(file_name)
# reproject to epsg:2180
df = df.to_crs(epsg=2180)

print(df.head())

# make a buffer of feature with name Świeradów Zdrój

feature = df[df['JPT_NAZWA_'] == 'Świeradów-Zdrój']

print(feature.head())

# make a buffer of feature with name Świeradów Zdrój
buffer = feature.buffer(200)
