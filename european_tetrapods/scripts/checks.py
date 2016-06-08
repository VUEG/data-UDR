import glob
import logging
import numpy as np
import numpy.ma as ma
import pprint
import rasterio
import sys
import yaml


# Check that the extent of all translated files matches and that all files
# have non-missing and non-zero values

# Helper function to check equality of list items
def add_item(dic, key, value):
    if key in dic.keys():
        dic[key] = dic[key] + [value]
    else:
        dic[key] = [value]
    return dic

# Get the inputs and outputs
input_rasters = snakemake.input
log_file = snakemake.log[0]

logger = logging.getLogger("check_data")
logger.setLevel(logging.DEBUG)
logFormatter = logging.Formatter("%(asctime)s [%(name)-12s] " +
                                 "[%(levelname)-5.5s] %(message)s",
                                 datefmt='%a, %d %b %Y %H:%M:%S')

fileHandler = logging.FileHandler(log_file, mode='w')
consoleHandler = logging.StreamHandler()
fileHandler.setFormatter(logFormatter)
logger.addHandler(fileHandler)
logger.addHandler(consoleHandler)


if len(input_rasters) == 0:
    print("No rasters found in the input folder")
    sys.exit(0)

# Make lists to holds the bound values
minx = {}
miny = {}
maxx = {}
maxy = {}
# Rasters with all values missing
empty_rasters = []
# Rasters with all values zero
zero_rasters = []

n_rasters = len(input_rasters)

for i, raster in enumerate(input_rasters):
    print('Checking raster {0}/{1}'.format(i+1, n_rasters), end="\r")
    with rasterio.open(raster) as src:
        minx = add_item(minx, src.bounds[0], raster)
        miny = add_item(miny, src.bounds[1], raster)
        maxx = add_item(maxx, src.bounds[2], raster)
        maxy = add_item(maxy, src.bounds[3], raster)

        src_data = src.read(1, masked=True)

        # Are all cells NoData?
        if ma.getmask(src_data).all():
            empty_rasters.append(raster)
        # Are all cells zero?
        if src_data.min() == 0.0 and src_data.max() == 0.0:
            zero_rasters.append(raster)
print()

print("{0} rasters found".format(len(input_rasters)))

if len(minx.keys()) > 1:
    logger.warning("Multiple minx found")
    pprint.pprint(minx)
else:
    logger.info("Only one minx found: {0}".format(list(minx)[0]))
if len(miny.keys()) > 1:
    logger.warning("Multiple miny found")
    pprint.pprint(miny)
else:
    logger.info("Only one miny found: {0}".format(list(miny)[0]))
if len(maxx.keys()) > 1:
    logger.warning("Multiple maxx found")
    pprint.pprint(maxx)
else:
    logger.info("Only one maxx found: {0}".format(list(maxx)[0]))
if len(maxy.keys()) > 1:
    logger.warning("Multiple maxy found")
    pprint.pprint(maxy)
else:
    logger.info("Only one maxy found: {0}".format(list(maxy)[0]))

if len(empty_rasters) > 0:
    logger.warning("Empty rasters found")
    for raster in empty_rasters:
        print("  " + raster)

if len(zero_rasters) > 0:
    logger.warning("Zero rasters found:")
    for raster in zero_rasters:
        logger.warning(" " + raster)

with open(snakemake.output[0], 'w') as outfile:
    outfile.write(yaml.dump(empty_rasters + zero_rasters,
                            default_flow_style=False))
