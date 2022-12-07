# -*- coding: utf-8 -*-
from dask.distributed import LocalCluster,Client
import stackstac
import dask.array as da
import pystac_client
import geopandas as gpd
import planetary_computer as pc
import xarray as xr

import os

from PIL import Image, ImageDraw, ImageFont
from xarray.plot.utils import _rescale_imshow_rgb
import numpy as np
from typing import TYPE_CHECKING, BinaryIO, cast,Literal
import ntpath

from IPython.display import HTML, display
import folium
import folium.plugins
from branca.element import Figure
import shapely.geometry
import matplotlib.pyplot as plt
import rioxarray
import pandas as pd
import rich.table

from dask.utils import ensure_dict, format_bytes

pc.settings.set_subscription_key('5f76969537404818891dee8872288a72')

cluster = LocalCluster(n_workers=5,
                       threads_per_worker=2,
                       dashboard_address=8787,
                       memory_limit='4GB')

client = Client(cluster)
#display(client)

wk = client.scheduler_info()["workers"]

text="Workers= " + str(len(wk))
memory = [w["memory_limit"] for w in wk.values()]
cores = sum(w["nthreads"] for w in wk.values())
text += ", Cores=" + str(cores)
if all(memory):
    text += ", Memory=" + format_bytes(sum(memory))
print(text)

# Define general variables

setting = {'spatial_resolution' : 100,
            'cloud_cover': 1}

glacier = "argentiere"
full_name= "Argentiere"

def convert_bounds(bbox, invert_y=False):
    """
    Helper method for changing bounding box representation to leaflet notation
    ``(lon1, lat1, lon2, lat2) -> ((lat1, lon1), (lat2, lon2))``
    """
    x1, y1, x2, y2 = bbox
    if invert_y:
        y1, y2 = y2, y1
    return ((y1, x1), (y2, x2))


# change directory 
os.chdir("C:/Users/kneibm/Documents/CAIRN/Remote_sensing/Avalanche_mapping/Planetary/earthgif/")

#filename = "geojson file path"
filename = glacier +".geojson"
# read in AOI as a GeoDataFrame
aoi = gpd.read_file(filename)

bbox =aoi.unary_union.bounds

ImSatSTAC = pystac_client.Client.open('https://planetarycomputer.microsoft.com/api/stac/v1')

#for collection in LandsatSTAC.get_collections():
    #print(collection)

search = (
    ImSatSTAC
    .search(
        bbox=bbox,
        datetime = "2019-01-01/2019-12-30", 
        collections = ["sentinel-1-rtc"],
    )
)
items = pc.sign(search)
print(str(len(items))+ ' scenes found')

stack = stackstac.stack(items, bounds_latlon=bbox, epsg = 32632,resolution=10)
stack

#stackstac.show(stack, center=None, zoom=None, range=None, cmap=None, checkerboard=True, interpolation='linear')


stack['time'].values
stack['s1:orbit_source'].values
stack['s1:processing_level'].values
stack['sar:product_type'].values
stack['sat:relative_orbit'].values
stack['sat:orbit_state'].values
stack['sar:instrument_mode'].values
stack.dims
stack.attrs

item = items[0]
table = rich.table.Table("key", "value")
for k, v in sorted(item.properties.items()):
    table.add_row(k, str(v))
table

VV = stack.sel(band="vv")
VV_Desc = VV.sel(s1:orbit_state="descending")


imgSelection = stack.sel(band="vv")[0].compute()

vv = stack.sel(band="vv")[0].compute()
#vv.plot.hist(bins=30);

# values need to be rescaled to be visualized.
def db_scale(x):
    return 10 * np.log10(x)

db_scale(vv).plot.hist(bins=50)

img = (
    db_scale(vv)
    .coarsen(x=1, y=1, boundary="trim")
    .max()
    .plot.imshow(cmap="bone", size=8,  add_colorbar=False)
)
img.axes.set_axis_off();


vh = stack.sel(band="vh")[0].compute()
img = (
    db_scale(vh)
    .coarsen(x=1, y=1, boundary="trim")
    .max()
    .plot.imshow(cmap="bone", size=8, add_colorbar=False)
)
img.axes.set_axis_off();

# export as geotiff
db_scale(vh).rio.to_raster(raster_path="test.tif")



