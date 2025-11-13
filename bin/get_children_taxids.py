#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import requests
from tqdm import tqdm
import sys
import argparse
import xml.etree.ElementTree as ET
import logging

from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_delay,
    wait_exponential,
    before_sleep_log,
)

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

# Modern NCBI API
NCBI_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
NCBI_API_HEADERS = {"accept": "application/json", "content-type": "application/json"}

# E-UTILITIES OL API
ESEARCH_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESEARCH_RETMAX = 1000000000  # max retmax that worked

TAXID_TO_NAME_OUTFILE_SUFFIX = ".taxids2names.csv"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute general statistics from count data for each sample"
    )
    parser.add_argument("--family", type=str, required=True, help="Family name")
    parser.add_argument(
        "--ncbi-api-key", dest="ncbi_api_key", type=str, help="NCBI API key"
    )
    return parser.parse_args()


@retry(
    retry=retry_if_exception_type(requests.exceptions.HTTPError),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_taxonomy(taxid: str | int, ncbi_api_key: str | None):
    taxons = [str(taxid)]
    data = {"taxons": taxons}
    headers = dict(NCBI_API_HEADERS)
    if ncbi_api_key is not None:
        headers["api_key"] = ncbi_api_key
    response = requests.post(NCBI_API_URL, headers=headers, json=data)
    response.raise_for_status()
    return response.json()


@retry(
    retry=retry_if_exception_type(requests.exceptions.HTTPError),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_esearch_query(query: str, database: str, ncbi_api_key: str | None):
    """
    Query NCBI's db with Esearch API
    """
    params = dict(db=database, term=query, retmax=ESEARCH_RETMAX)
    if ncbi_api_key is not None:
        params["api_key"] = ncbi_api_key
    response = requests.get(ESEARCH_BASE_URL, params=params)
    response.raise_for_status()
    return response.text


def parse_ids_from_xml(xml_string: str):
    """
    Parse XML string and get text in <Id>...</Id> blocks
    :param xml_string:
    :return: list of Ids
    """
    # parsing XML returned by API
    root = ET.fromstring(xml_string)
    # species taxon IDs are contained in <Id>...</Id> blocks
    return [id_element.text for id_element in root.findall(".//Id")]


def get_taxon_metadata(taxid: str, ncbi_api_key: str) -> dict:
    result = send_request_to_ncbi_taxonomy(taxid, ncbi_api_key)
    if len(result["taxonomy_nodes"]) > 1:
        raise ValueError(f"Multiple taxids for family {family}")
    return result["taxonomy_nodes"][0]


def get_family_taxid(family: str, ncbi_api_key: str):
    metadata = get_taxon_metadata(family, ncbi_api_key)
    if "taxonomy" not in metadata:
        logger.info(f"Could not find taxonomy results for family {family}")
        if "errors" in metadata:
            for error in metadata["errors"]:
                logger.error(f"Error: {error['reason']}\n")
                sys.exit(100)
    return int(metadata["taxonomy"]["tax_id"])


def get_children_taxids(taxid: str, ncbi_api_key: str | None) -> list[str]:
    """
    Get list of children taxonomy IDs given a family taxonomy ID
    :param taxid:
    :return: list of SRA experiment IDs
    """
    xml_string = send_esearch_query(
        query=f"txid{taxid}[Subtree]",
        # query=f"txid{taxid}[Subtree] AND species[Rank]",
        database="taxonomy",
        ncbi_api_key=ncbi_api_key,
    )
    return parse_ids_from_xml(xml_string)


def get_taxon_name(taxid: str):
    metadata = get_taxon_metadata(taxid)
    return metadata.get("taxonomy", {}).get("organism_name", None)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()
    family = args.family

    family_taxid = get_family_taxid(family, args.ncbi_api_key)
    logger.info(f"Family taxid: {family_taxid}")

    logger.info(f"Getting children taxids for family {family}")
    children_taxids = get_children_taxids(family_taxid, args.ncbi_api_key)

    # converting all taxids to int for uniformity
    # adding family taxid to children taxids (in case)
    # keeping unique taxids (in case)
    children_taxids = list(
        set([int(taxid) for taxid in children_taxids + [family_taxid]])
    )
    logger.info(f"Obtained {len(children_taxids)} children taxids\n")

    taxids2names = {taxid: get_taxon_name(taxid) for taxid in tqdm(children_taxids)}

    taxid2name_outfile = f"{family}{TAXID_TO_NAME_OUTFILE_SUFFIX}"
    with open(taxid2name_outfile, "w") as fout:
        for taxid, name in taxids2names.items():
            fout.write(f"{taxid},{name}\n")

    logger.info("Done")
