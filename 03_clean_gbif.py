import argparse
import io
import os
from pathlib import Path
from typing import Dict, List, Union

import requests
import zipfile
import pandas as pd

DOI = "10.15468/dl.pan5fn"


def download_gbif_doi(doi: str, output_path: Path) -> None:
    response = requests.get(f"https://api.gbif.org/v1/occurrence/download/{doi}")
    download_key = response.json()['key']
    response = requests.get(f"https://api.gbif.org/v1/occurrence/download/request/{download_key}")

    with zipfile.ZipFile(io.BytesIO(response.content)) as zip_file:
        file_list = zip_file.namelist()
        csv_files = [f for f in file_list if f.endswith('.csv')]
        csv_filename = csv_files[0]
        with zip_file.open(csv_filename) as csv_file:
            with open(output_path, 'wb') as output_file:
                output_file.write(csv_file.read())


if not Path("data", "raw_occurrence.csv").exists():
    download_gbif_doi(DOI, Path("data", "raw_occurrence.csv"))

occ_df = pd.read_csv(Path("data", "raw_occurrence.csv"), delimiter='\t')
key_df = pd.read_csv(Path("data", "keys.csv"))

column_mapping = {
    'species' : 'species',
    'decimalLatitude': 'latitude',
    'decimalLongitude': 'longitude'
}

os.mkdir(Path("data", "species_occurrence"))

for (i, row) in key_df.iterrows():
    this_species = occ_df[occ_df.taxonKey == row.key]    
    this_species = this_species.rename(columns=column_mapping)
    this_species = this_species[['latitude', 'longitude']]
    this_species = this_species.dropna(subset=['latitude', 'longitude'])
    this_species.to_csv(Path("data", "species_occurrence", f"{row.taxon}.csv"))

