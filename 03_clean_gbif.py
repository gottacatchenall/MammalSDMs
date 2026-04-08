import argparse
import io
import os
from pathlib import Path
from typing import Dict, List, Union

import pandas as pd


occ_df = pd.read_csv(Path("data", "raw_occurrence.csv"), delimiter='\t')
key_df = pd.read_csv(Path("data", "keys.csv"))


column_mapping = {
    'species' : 'species',
    'decimalLatitude': 'latitude',
    'decimalLongitude': 'longitude'
}

for (i, row) in key_df.iterrows():
    this_species = occ_df[occ_df.taxonKey == row.key]    
    this_species = this_species.rename(columns=column_mapping)
    this_species = this_species[['latitude', 'longitude']]
    this_species = this_species.dropna(subset=['latitude', 'longitude'])
    this_species.to_csv(Path("data", "species_occurrence", f"{row.taxon}.csv"))

