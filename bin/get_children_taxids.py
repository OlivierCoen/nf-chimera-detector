#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
Find all children (here only species and eerything below) taxids, given a taxid
"""

import argparse
import logging
import sys
import json
import xml.etree.ElementTree as ET

import requests
from tenacity import (
    before_sleep_log,
    retry,
    retry_if_exception_type,
    stop_after_delay,
    wait_exponential,
)
from tqdm import tqdm

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

# Modern NCBI API
NCBI_API_DATASET_REPORT_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/dataset_report"
NCBI_API_HEADERS = {"accept": "application/json", "content-type": "application/json"}
NCBI_DATASET_REPORT_TAXON_CHUNKSIZE = 1000

# E-UTILITIES OL API
ESEARCH_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESEARCH_RETMAX = 1000000000  # max retmax that worked

CHILDREN_TAXON_METADATA_OUTFILE_SUFFIX = ".children_taxons_metadata.json"
TAXID_TO_NAME_OUTFILE_SUFFIX = ".taxids2names.csv"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get all descendant NCBI taxon IDs from a specific taxon ID"
    )
    parser.add_argument("--taxon", type=str, required=True, help="Taxon name")
    parser.add_argument("--keep-below", dest="upper_node_type_allowed", type=str, required=True, help="Keep all nodes below this type of node")
    return parser.parse_args()

class RateLimitException(Exception):
    pass

@retry(
    retry=retry_if_exception_type((
        RateLimitException,
        requests.exceptions.HTTPError
    )),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def get_metadata_for_taxons(taxons: list[str]):
    chunked_taxons_list = [
        taxons[i:i + NCBI_DATASET_REPORT_TAXON_CHUNKSIZE] 
        for i in range(0, len(taxons), NCBI_DATASET_REPORT_TAXON_CHUNKSIZE)
    ]
    for chunked_taxons in chunked_taxons_list:
        response = requests.post(
            NCBI_API_DATASET_REPORT_URL,
            headers=NCBI_API_HEADERS,
            json={'taxons': chunked_taxons}
        )
        if response.status_code == 429:
            raise RateLimitException("Rate limit exceeded")
        response.raise_for_status()
        yield response.json()


@retry(
    retry=retry_if_exception_type(requests.exceptions.HTTPError),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_esearch_query(query: str, database: str):
    """
    Query NCBI's db with Esearch API
    """
    params = dict(db=database, term=query, retmax=ESEARCH_RETMAX)
    response = requests.get(ESEARCH_BASE_URL, params=params)
    response.raise_for_status()
    return response.text


def parse_ids_from_xml(xml_string: str) -> list[str]:
    """
    Parse XML string and get text in <Id>...</Id> blocks
    :param xml_string:
    :return: list of Ids
    """
    # parsing XML returned by API
    root = ET.fromstring(xml_string)
    # species taxon IDs are contained in <Id>...</Id> blocks
    return [
        id_element.text
        for id_element in root.findall(".//Id")
        if id_element.text is not None
    ]


def get_parent_taxid(taxon: str) -> int:
    result = next(get_metadata_for_taxons([taxon]))
    if len(result["reports"]) > 1:
        raise ValueError(f"Multiple taxids for taxon {taxon}")
    metadata = result["reports"][0]
    if "taxonomy" not in metadata:
        logger.info(f"Could not find taxonomy results for taxon {taxon}")
        if "errors" in metadata:
            for error in metadata["errors"]:
                logger.error(f"Error: {error['reason']}\n")
                sys.exit(100)
    return int(metadata["taxonomy"]["tax_id"])


def get_all_children_taxids(taxid: int) -> list[str]:
    """
    Get list of all children taxonomy IDs given a taxon taxonomy ID
    :param taxid:
    :return: list of children IDs
    """
    xml_string = send_esearch_query(
        query=f"txid{taxid}[Subtree]",
        # query=f"txid{taxid}[Subtree] AND species[Rank]",
        database="taxonomy",
    )
    return parse_ids_from_xml(xml_string)


def is_valid_taxon(taxonomy_report: dict, upper_node_type_allowed: str) -> bool:
    # if one of the parents node is at the species level, we keep it
    if taxonomy_report.get('classification', {}).get(upper_node_type_allowed, {}).get('id') is not None:
        return True
    # handle exceptions
    # if rank is None ('no rank'), it most generally (but not always) means that it's under species, or at the same level
    if taxonomy_report.get('rank') is None:
        return True
    return False
    

def filter_children_taxids(taxids: list[int], upper_node_type_allowed: str) -> list[dict]:
    valid_taxons = []
    for chunk_result in get_metadata_for_taxons([str(taxid) for taxid in taxids]):
        for report in chunk_result['reports']:
            if 'taxonomy' not in report:
                continue
            taxonomy_report = report['taxonomy']
            if is_valid_taxon(taxonomy_report, upper_node_type_allowed):
                taxon = {
                    'taxid': taxonomy_report['tax_id'],
                    'name': taxonomy_report.get('current_scientific_name', {}).get('name', 'Unknown'),
                    'rank': taxonomy_report.get('rank', "RANK_UNKNOWN"),
                    'classification': taxonomy_report.get('classification', {})
                }
                valid_taxons.append(taxon)
    return valid_taxons


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()
    parent_taxon = args.taxon

    parent_taxid = get_parent_taxid(parent_taxon)
    logger.info(f"Parent taxid: {parent_taxid}")

    logger.info(f"Getting children taxids for family {parent_taxon}")
    children_taxids = get_all_children_taxids(parent_taxid)

    # converting all taxids to int for uniformity
    # adding family taxid to children taxids (in case)
    # keeping unique taxids (in case)
    children_taxids = sorted({int(taxid) for taxid in children_taxids + [parent_taxid]})
    logger.info(f"Obtained {len(children_taxids)} children taxids\n")

    logger.info(f"Filtering out all nodes above {args.upper_node_type_allowed} level")
    filtered_taxons = filter_children_taxids(children_taxids, args.upper_node_type_allowed)
    logger.info(f"Kept {len(filtered_taxons)} children taxons of species rank or below")

    children_taxon_metadata_outfile = f"{parent_taxon}{CHILDREN_TAXON_METADATA_OUTFILE_SUFFIX}"
    with open(children_taxon_metadata_outfile, "w") as fout:
        json.dump(filtered_taxons, fout)

    taxid2name_outfile = f"{parent_taxon}{TAXID_TO_NAME_OUTFILE_SUFFIX}"
    with open(taxid2name_outfile, "w") as fout:
        for taxon in filtered_taxons:
            fout.write(f"{taxon['taxid']},{taxon['name']}\n")

    logger.info("Done")
