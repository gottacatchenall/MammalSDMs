import argparse
import io
import json
import os
import requests
import subprocess
import zipfile
from pathlib import Path
from typing import Dict, List, Union

import pandas as pd
import pygbif


class GBIFQueryBuilder:
    """Build and execute GBIF occurrence queries."""
    GBIF_API_URL = "https://api.gbif.org/v1/occurrence/download/request"

    def __init__(self, template_path: Path, credentials: Dict[str, str], key_path: Path = Path("data", "keys.csv")):
        """Initialize query builder.

        Args:
            template_path: Path to GBIF query template JSON
            credentials: Dict with 'username' and 'password' keys
        """
        assert template_path.exists(), f"Template not found: {template_path}"
        assert credentials.get('username') and credentials.get('password'), \
            "Credentials must include 'username' and 'password'"

        with open(template_path, 'r') as f:
            self.template = json.load(f)

        self.credentials = credentials
        self.key_path = key_path

    def _get_taxon_keys(self, taxa: List[str]) -> Dict[str, int]:
        taxon_keys = {}
        for taxon in taxa:
            response = pygbif.species.name_backbone(
                scientificName=taxon,
                taxonRank="species"
            )
            key = response.get('usage', {}).get('key')
            print(f"{taxon}: {key}")
            assert key is not None, f"Could not find GBIF key for: {taxon}"
            taxon_keys[taxon] = int(key)
        return taxon_keys

    def add_taxon_filter(self, taxa: List[str]) -> None:
        taxon_keys = self._get_taxon_keys(taxa)
        pd.DataFrame([[a,b] for a,b in taxon_keys.items()], columns = ["taxon", "key"]).to_csv(Path("data", "keys.csv"))
        self.template['predicate']['predicates'].append({
            "type": "in",
            "key": "TAXON_KEY",
            "values": list(taxon_keys.values())
        })

    def set_creator(self, username: str) -> None:
        self.template['creator'] = username

    def get_query(self) -> Dict:
        return self.template

    def execute(self, debug: bool = True) -> str:
        if debug:
            print("Request payload:")
            print(json.dumps(self.template, indent=2))
            print()

        cmd = [
            "curl",
            "--include",
            "--user", f"{self.credentials['username']}:{self.credentials['password']}",
            "--header", "Content-Type: application/json",
            "--data", json.dumps(self.template),
            self.GBIF_API_URL
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"GBIF query failed: {result.stderr}"
        return result.stdout


class EnvironmentManager:
    @staticmethod
    def load_env_file(dotenv_path: Path = Path(".env")) -> None:
        if not dotenv_path.exists():
            return

        with open(dotenv_path, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    key, value = line.strip().split('=', 1)
                    os.environ[key] = value

    @staticmethod
    def get_gbif_credentials() -> Dict[str, str]:
        username = os.environ.get("GBIF_USERNAME")
        password = os.environ.get("GBIF_PASSWORD")
        assert username and password, \
            "GBIF_USERNAME and GBIF_PASSWORD must be set in environment"

        return {"username": username, "password": password}


class GBIFOccurrenceDownloader:
    def __init__(
        self,
        template_path: Path = Path("data") / "raw" / "Ticks" / "gbif_query_template.json",
        env_path: Path = Path(".env")
    ):
        self.template_path = Path(template_path)

        EnvironmentManager.load_env_file(env_path)
        self.credentials = EnvironmentManager.get_gbif_credentials()

    def query(self, taxa: List[str], debug: bool = True) -> str:
        # Build query
        builder = GBIFQueryBuilder(self.template_path, self.credentials)
        builder.add_taxon_filter(taxa)

        # Set creator
        builder.set_creator(self.credentials['username'])

        # Execute
        return builder.execute(debug=debug)


def main():
    parser = argparse.ArgumentParser(
        description="Download GBIF occurrence data for specified taxa.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--taxa_list",
        nargs="+",
        default=Path("data", "species_list.txt"),
        help="List of scientific names to query"
    )
    parser.add_argument(
        "--template",
        nargs="+",
        default=Path("data", "gbif_query_template.json"),
        help="Template JSON for GBIF query"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=True,
        help="Print request payload before executing query (default: True)"
    )

    args = parser.parse_args()
    try:
        downloader = GBIFOccurrenceDownloader(template_path=args.template)
        with open(args.taxa_list) as f:
            taxa = f.read().splitlines()
            response = downloader.query(taxa, debug=args.debug)
            print(response)
    except AssertionError as e:
        print(f"Error: {e}", file=__import__('sys').stderr)
        exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=__import__('sys').stderr)
        exit(1)


def download_gbif_doi(doi: str, output_path: Path) -> None:
    # Get download key from DOI
    response = requests.get(f"https://api.gbif.org/v1/occurrence/download/{doi}")
    download_key = response.json()['key']

    # Download the zip file
    response = requests.get(f"https://api.gbif.org/v1/occurrence/download/request/{download_key}")

    # Parse the zip file from response content
    with zipfile.ZipFile(io.BytesIO(response.content)) as zip_file:
        file_list = zip_file.namelist()
        csv_files = [f for f in file_list if f.endswith('.csv')]

        if not csv_files:
            raise ValueError("No CSV files found in the zip archive")

        # Extract the first CSV file
        csv_filename = csv_files[0]

        # Extract and save to output_path
        with zip_file.open(csv_filename) as csv_file:
            with open(output_path, 'wb') as output_file:
                output_file.write(csv_file.read())

    new_df = clean_gbif_occurrences(pd.read_csv(output_path))
    new_df.to_csv(output_path, index=False)


def clean_gbif_occurrences(df: pd.DataFrame) -> pd.DataFrame:
    cleaned_df = df.copy()
    column_mapping = {
        'species' : 'species',
        'decimalLatitude': 'latitude',
        'decimalLongitude': 'longitude'
    }
    cleaned_df = cleaned_df.rename(columns=column_mapping)

    # Select only the columns we need
    cleaned_df = cleaned_df[['species', 'latitude', 'longitude']]

    # Drop rows with missing coordinates
    cleaned_df = cleaned_df.dropna(subset=['latitude', 'longitude'])
    return cleaned_df


if __name__ == "__main__":
    main()