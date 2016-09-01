### Folder structure

```
├── amphibians          <- Processed data for amphibians.
│
├── birds               <- Processed data for birds.
│
├── mammals             <- Processed data for mammals.
│
├── reptiles            <- Original data for reptiles.
│
├── data_manifest.yml   <- A YAML file listing all the processed and included
│                          files. This files can serve as a download manifest
│                          file for any project using the dataset.
│
└── spp_codes.csv       <- A CSV file containing a lookup table linking the
                           numeric codes in the AIG names into real species
                           names.
```

### Raster naming convention

All files have the species' scientific name in them.

### Cell values

For each species `j` and cell `i`, the value `r_ij` is in range [0, 100] and
corresponds to the estimated percentage of primary and marginal habitat type 
for the given species.
