## Distributions of European tetrapods

### Description

This dataset include the distributions of 819 species of European tetrapods (amphibians, birds, mammals, and reptiles). For the four groups, the EOOs (extent of occurrence) were then refined using habitat preferences for all species, obtained from expert opinion and published literature. The resolution of the data is 1 square kilometer.

In the original European tetrapods dataset, there are 819 species. For the version here, the original data has been subset spatially to include EU-26 countries (AT, BE, BG, CZ, DE, DK, ES, EL, EE, FR, FI, IT, HU, IE, NL, LU, LI, LT, LV, PL, SE, RO, PT, SK, SI, UK). After the spatial subsetting, this dataset contains 763 species. The number of species in each species group is:

| Group      | N       |
|------------|---------|
| Amphibians | 83      |
| Birds      | 404     |
| Mammals    | 164     |
| Reptiles   | 112     |
| **Total**  | **763** |

See also https://github.com/VUEG/data-UDR for a development version.

### Folder structure and processing

```
├── org                 <- Original data received from Luigi Maiorano       
│                          <luigi.maiorano@uniroma1.it> on 2016-08-03 via  
│                          Dropbox. Copied manually from here since couldn't get  
│                          snakemake's Dropbox remote working.  
└── european_tetrapods  <- Processed and filtered files
```

The files in sub-directory `org` are as they have been received from Luigi Maiorano. Folder `european_tetrapods` contains processed files. The processing steps are the following:

1. Check data for empty layers
2. Translate numeric codes into proper species names
3. Translate all rasters from ArcInfo grids (AIGs) to GeoTIFFs.

See [Snakemake file](https://github.com/VUEG/data-UDR/blob/master/european_tetrapods/Snakefile) for details.

### Metadata

#### Cell values

For each species `j` and cell `i`, the value `r_ij` is in range [0, 100] and corresponds to the estimated pSanekmakeercentage of *primary* and *margnial* habitat for the given species. Original data can separate between the two habitat types, but here they have been pooled together.

#### Geospatial

Following metadata applies to all individual species rasters in the dataset.

```
"res": [
  1000.0,
  1000.0
],
"bounds": [
  1999999.999999999,
  1000000.0000000009,
  6525999.999999999,
  5410000.000000001
],
"dtype": "uint8",
"driver": "GTiff",
"transform": [
  1000.0,
  0.0,
  1999999.999999999,
  0.0,
  -1000.0,
  5410000.000000001
],
"lnglat": [
  4262999.999999999,
  3205000.000000001
],
"height": 4410,
"width": 4526,
"shape": [
  4410,
  4526
],
"blockxsize": 4526,
"tiled": false,
"blockysize": 1,
  "nodata": 255.0

```

**NOTE:** CRS has not explicitly defined for the rasters, but it is ETRS-LAEA
(EPSG:3035).

### License

All rights reserved. For permission to use the data, contact Luigi Maiorano (<luigi.maiorano@uniroma1.it>) at the University of Rome.

### Contributors

+ Luigi Maiorano (<luigi.maiorano@uniroma1.it>)
+ Joona Lehtomäki (<joona.lehtomaki@gmail.com>)

### References

+ Maiorano, L., Amori, G., Capula, M., Falcucci, A., Masi, M., Montemaggiori, A., … Guisan, A. (2013). Threats from Climate Change to Terrestrial Vertebrate Hotspots in Europe. PLoS ONE, 8(9), 1–14. http://doi.org/10.1371/journal.pone.0074989
+ Thuiller, W., Maiorano, L., Mazel, F., Guilhaumon, F., Ficetola, G. F., Lavergne, S., … Mouillot, D. (2015). Conserving the functional and phylogenetic trees of life of European tetrapods. Philosophical Transactions of the Royal Society of London. Series B, Biological Sciences, 370(1662), 20140005. http://doi.org/10.1098/rstb.2014.0005
