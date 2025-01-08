import geopandas as gpd
import requests
import pathlib

# read borders
border = gpd.read_file("projekt_2_granice/strefa wielkomiejska_bufor_200.shp").to_crs("EPSG:4326")

# check if osm_data folder exists
osm_data_folder = "osm_data"
pathlib.Path(osm_data_folder).mkdir(parents=True, exist_ok=True)

# check if lodzkie-latest.osm.pbf exists
if not pathlib.Path("lodzkie-latest.osm.pbf").exists():
    # download the file
    url = "https://download.geofabrik.de/europe/poland/lodzkie-latest.osm.pbf"
    r = requests.get(url, allow_redirects=True)
    with open("lodzkie-latest.osm.pbf", "wb") as f:
        f.write(r.content)

# read the data
for idx, layer in gpd.list_layers("lodzkie-latest.osm.pbf").iterrows():
    name = layer["name"]
    print(name)
    df = gpd.read_file("lodzkie-latest.osm.pbf", layer=name)
    print(df.head())
    # cut to the border and save as name.gpkg
    cut = df.clip(border).to_crs("EPSG:2180")
    cut.to_file(f"{osm_data_folder}/{name}.gpkg")

