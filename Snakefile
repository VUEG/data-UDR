# Rule to checkt the extents of tranlsated bird files.
rule check_extents:
    input:
        "data/birds"
    output:
        "data/birds"
    log:
        "log/extents.log"
    message: "Checking bird raster extents..."
    script:
        "scripts/check_extents.py"

# Rule to extract the data from from the 7z archive
rule extract:
    input:
        "data/originals/birds.7z"
    output:
        "data"
    log:
        "log/extract.log"
    message: "Extracting AIGs from {input}."
    shell:
        # Extra quotes needed because of the whitespace in the name.
        # NOTE: This assumes that 7zip is around.
        # NOTE: Only subfolder birds/output and ReadMe.xlsx are actually extracted.
        "7za x '{input}' birds/output birds/ReadMe.xlsx -o{output} -r -aoa -mmt{threads} >& {log}"

# Rule to convert the Excel file containing the species names into CSV.
rule spconvert:
    input:
        "data/birds/ReadMe.xlsx"
    output:
        "data/birds/farmland_birds_sp.csv"
    log:
        "log/spconvert.log"
    message: "converting {input} to CSV."
    script:
        "scripts/spconvert.py"

# Rule to translate a filtered set of bird species distributions from AIG to
# GeoTIFFs. File names will be replaced with the species name.
rule sptranslate:
    input:
        input_dir='data/birds/output',
        sp_data='data/birds/farmland_birds_sp.csv'
    output:
        "data/birds"
    log:
        "log/sptranslate.log"
    message: "Translating files to GeoTIFF"
    script:
        "scripts/sptranslate.py"
