### Folder structure

```
├── amph                <- Original data for amphibians.
│
├── bird                <- Original data for birds.
│
├── mam                 <- Original data for mammals.
│
├── rept                <- Original data for reptiles.
│
├── bird                <- Original data for birds.
│
├── Codes.xlsx          <- Excel worksheet lookup table linking the numeric
│                          codes in the AIG names into real species names.
│
└── skpped_rasters.yml  <- A YAML file listing all the empty rasters (these
                           are not translated over to european_tetrapods).
```

### Raster naming convention

Rasters here are ArcInfo grids (AIGs). For each species, there are two rasters:

```
m[a|b|m|r][NO]_[2|12]
```

Here, `m` doesn't stand for anything in particular, `a|b|m|r` stands for the
species group (`a`mphibians, `b`irds, `m`ammals, `r`eptiles), `NO` for a
species-specific numeric code, and `2|12` for which habitat types are included
(`12` = primary and marginal, `2` = marginal).

### Cell values

For each species `j` and cell `i`, the value `r_ij` is in range [0, 100] and corresponds to the estimated percentage of a given habitat type for the given species.
