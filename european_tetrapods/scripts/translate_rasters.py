import errno
import logging
import os
import subprocess
import sys
import re
import yaml
from pandas import read_csv


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


# Get the inputs and outputs
input_rasters = snakemake.input['datasets']
spp_names_file = snakemake.input['spp_names_file']
skip_rasters = snakemake.input.skip_rasters
output_dirs = snakemake.output
log_file = snakemake.log[0]

logger = logging.getLogger("translate_rasters")
logger.setLevel(logging.DEBUG)
logFormatter = logging.Formatter("%(asctime)s [%(name)-12s] " +
                                 "[%(levelname)-5.5s] %(message)s",
                                 datefmt='%Y-%m-%d %H:%M:%S')

fileHandler = logging.FileHandler(log_file, mode='w')
consoleHandler = logging.StreamHandler()
fileHandler.setFormatter(logFormatter)
logger.addHandler(fileHandler)
logger.addHandler(consoleHandler)

SKIP_RASTERS = yaml.safe_load(open(skip_rasters, 'r'))


# Load the additional CSV data
spp_data = read_csv(spp_names_file)

n_files = len(input_rasters)
counter = 1
# Create a data manifest dictionary
data_manifest = {'provider': 'udr',
                 'collections': {'amphibians': [],
                                 'birds': [],
                                 'mammals': [],
                                 'reptiles': []}}

for aig in input_rasters:
    # Skip raster if it's in skip_rasters list
    if aig in SKIP_RASTERS:
        logger.warning('[{0}/{1}] Skipping an empty'.format(counter, n_files) +
                       ' or zero raster {}'.format(aig))
        next
    # Figure out the species. First, get rid of 'hdr.adf' in the AIG name
    spp_path = aig.replace('/hdr.adf', '')
    # Split the path
    path_components = spp_path.split('/')
    # Then, get the spp code (last element of the path, so use pop)
    spp_code = path_components[-1]
    # Fore reason or another, spp_code is prefix with an additional "m". Get
    # rid of that
    spp_code = spp_code[1:]

    # Get the auxiliary data
    sci_name = spp_data.loc[spp_data.code == spp_code, 'species']
    taxon = spp_data.loc[spp_data.code == spp_code, 'class'].values[0].lower()
    # Finally, fix the output path
    output_path = [path for path in output_dirs if taxon in path][0]

    if len(sci_name) == 0:
        logger.warning('[{0}/{1}] No entry found'.format(counter, n_files) +
                       ' for AIG {0} (spp_code: {1})'.format(aig, spp_code))
    else:
        sci_name = sci_name.values[0]
        output_file = sci_name.lower().replace(' ', '_') + '.tif'
        # Get rid of problematic characters in the names
        output_file = output_file.replace('(', '')
        output_file = output_file.replace(')', '')
        output_file = os.path.join(output_path, output_file)
        # Create the folder if it does not exist
        if not os.path.exists(os.path.dirname(output_file)):
            mkdir_p(os.path.dirname(output_file))
        # Define potentially platform-specific, GDAL-related variables
        # constants
        cmd = ['gdal_translate', '-of', 'GTiff', '-co', 'COMPRESS=DEFLATE']
        cmd.append(aig)
        cmd.append(output_file)
        logger.info('[{0}/{1}] Translating {2} to {3}...'.format(counter,
                                                                 n_files,
                                                                 aig,
                                                                 output_file))
        subprocess.run(cmd)
        file_base = os.path.basename(output_file)
        data_manifest['collections'][taxon].append(file_base)
        counter += 1

data_manifest['collections']['amphibians'].sort()
data_manifest['collections']['birds'].sort()
data_manifest['collections']['mammals'].sort()
data_manifest['collections']['reptiles'].sort()

# Write the data manifest file
outfile = snakemake.output['data_manifest']
with open(outfile, 'w') as outfile:
    outfile.write(yaml.dump(data_manifest, default_flow_style=False))
logger.info('Created data manifest file {}'.format(outfile))
