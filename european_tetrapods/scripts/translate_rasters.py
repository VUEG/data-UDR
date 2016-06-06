import errno
import logging
import os
import subprocess
import sys
import re
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
output_dirs = snakemake.output
log_file = snakemake.log[0]

logger = logging.getLogger("translate_rasters")
logger.setLevel(logging.DEBUG)
logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")

fileHandler = logging.FileHandler(log_file)
fileHandler.setFormatter(logFormatter)
logger.addHandler(fileHandler)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

# Load the additional CSV data
spp_data =  read_csv(spp_names_file)

n_files = len(input_rasters)
counter = 1

for aig in input_rasters:
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
        logger.warning('No entry found for AIG {0} (spp_code: {1})'.format(aig, spp_code))
    else:
        sci_name = sci_name.values[0]
        output_file = sci_name.lower().replace(' ', '_') + '.tif'
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
        counter += 1
