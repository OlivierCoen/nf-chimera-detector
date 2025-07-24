#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import requests
import sys
import argparse
import logging

from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_delay,
    wait_exponential,
    before_sleep_log,
)

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

# Modern NCBI API
NCBI_GENOME_DATASET_REPORT_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report"
NCBI_GENOME_DATASET_REPORT_API_PARAMS = "filters.assembly_level=chromosome&filters.assembly_level=complete_genome&page_size=1000"

NCBI_API_TAXONOMY_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"

NCBI_API_HEADERS = {
    "accept": "application/json",
    "content-type": "application/json"
}

OUTFILE = "mean_assembly_length.txt"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get assembly statistics for genomes corresponding to children taxids. "
                    "If no genome assembly was found in children taxids, research is recursively performed on parent taxids."
    )
    parser.add_argument(
        "--family", type=str, required=True, help="Family name"
    )
    return parser.parse_args()


@retry(
    retry=retry_if_exception_type(requests.exceptions.HTTPError),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_genome_dataset_api(taxid: int):

    url = NCBI_GENOME_DATASET_REPORT_API_URL.format(taxid=taxid)
    url += f"?{NCBI_GENOME_DATASET_REPORT_API_PARAMS}"

    response = requests.get(url, headers=NCBI_API_HEADERS)
    # check status
    if not response.ok:
        raise RuntimeError(f"Error: {response.status_code}. {response.text}")
    return response.json()


@retry(
    retry=retry_if_exception_type(requests.exceptions.HTTPError),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_taxonomy(family: str):
    taxons = [family]
    data = {
        "taxons": taxons
    }
    response = requests.post(NCBI_API_TAXONOMY_URL, headers=NCBI_API_HEADERS, json=data)
    # check status
    if not response.ok:
        raise RuntimeError(f"Error: {response.status_code}. {response.text}")
    return response.json()


def get_family_taxid(family: str):

    result = send_request_to_ncbi_taxonomy(family)

    if len(result['taxonomy_nodes']) > 1:
        raise ValueError(f"Multiple taxids for family {family}")
    node = result['taxonomy_nodes'][0]

    if "taxonomy" not in node:
        logger.info(f"Could not find taxonomy results for family {family}")
        if "errors" in node:
            for error in node["errors"]:
                logger.error(f"Error: {error['reason']}\n")
                sys.exit(100)

    return node['taxonomy']['tax_id']


def get_parent_taxid(taxid: int):

    result = send_request_to_ncbi_taxonomy(taxid)

    if 'taxonomy_nodes' not in result:
        raise ValueError(f"Could not get taxonomy data for taxid {taxid}")
    node = result['taxonomy_nodes'][0]

    if "taxonomy" not in node:
        logger.info(f"Could not find taxonomy results for taxid {taxid}")
        if "errors" in node:
            for error in node["errors"]:
                logger.error(f"Error: {error['reason']}\n")
                sys.exit(100)

    if not node["taxonomy"].get('lineage'):
        logger.info(f"Could not find parent for taxid {taxid}")
        sys.exit(100)

    return node['taxonomy']['lineage'][-1]



def get_mean_assembly_length(taxid: int, attempt: int) -> float:

    result = send_request_to_ncbi_genome_dataset_api(taxid)

    assembly_lengths = []
    if result.get('reports') and isinstance(result['reports'], list):
        for report in result['reports']:
            if assembly_length := report.get('assembly_stats', {}).get('total_sequence_length'):
                try:
                    assembly_length = int(assembly_length)
                except ValueError as err:
                    logger.error(f"Could not parse assembly length {assembly_length}")
                    continue
                assembly_lengths.append(assembly_length)

    # if, for any reason, we could not get any assembly stat, we try with parent
    if not assembly_lengths:
        logger.warning(f'Could not get assembly stats for children for taxid {taxid}. Trying with parent.')
        # recursive call with parent taxids
        parent_taxid = get_parent_taxid(taxid)
        logger.info(f"Parent taxid: {parent_taxid}")
        return get_mean_assembly_length(parent_taxid)

    # returning mean of assembly lengths (as integer)
    return int(sum(assembly_lengths) / len(assembly_lengths))


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()
    family = args.family

    family_taxid = get_family_taxid(family)
    logger.info(f"Family taxid: {family_taxid}")

    logger.info(f"Getting assembly stats for children in taxid {family_taxid}")

    mean_assembly_length = get_mean_assembly_length(family_taxid)
    logger.info(f"Mean assembly length: {mean_assembly_length}")

    with open(OUTFILE, 'w') as fout:
        fout.write(str(mean_assembly_length))

    logger.info("Done")
