#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import json
import logging
import sys

import requests
from tenacity import (
    before_sleep_log,
    retry,
    retry_if_exception_type,
    stop_after_delay,
    wait_exponential,
)

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

# Modern NCBI API
NCBI_GENOME_DATASET_REPORT_API_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report"
)
NCBI_GENOME_DATASET_REPORT_API_PARAMS = "filters.assembly_level=chromosome&filters.assembly_level=complete_genome&page_size=1000"

NCBI_API_TAXONOMY_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"

NCBI_API_HEADERS = {"accept": "application/json", "content-type": "application/json"}


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
    parser.add_argument("--family", type=str, required=True, help="Family name")
    parser.add_argument(
        "--out", dest="outfile", type=str, required=True, help="Outfile name"
    )
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
def send_request_to_ncbi_genome_dataset_api(taxid: int, ncbi_api_key: str | None):
    url = NCBI_GENOME_DATASET_REPORT_API_URL.format(taxid=taxid)
    url += f"?{NCBI_GENOME_DATASET_REPORT_API_PARAMS}"
    headers = dict(NCBI_API_HEADERS)
    if ncbi_api_key is not None:
        headers["api_key"] = ncbi_api_key
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    return response.json()


@retry(
    retry=retry_if_exception_type(requests.exceptions.HTTPError),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_taxonomy(taxon: str | int, ncbi_api_key: str | None):
    taxons = [str(taxon)]
    data = {"taxons": taxons}
    headers = dict(NCBI_API_HEADERS)
    if ncbi_api_key is not None:
        headers["api_key"] = ncbi_api_key
    response = requests.post(NCBI_API_TAXONOMY_URL, headers=headers, json=data)
    response.raise_for_status()
    return response.json()


def get_family_taxid(family: str, ncbi_api_key: str | None):
    result = send_request_to_ncbi_taxonomy(family, ncbi_api_key)

    if len(result["taxonomy_nodes"]) > 1:
        raise ValueError(f"Multiple taxids for family {family}")
    node = result["taxonomy_nodes"][0]

    if "taxonomy" not in node:
        logger.info(f"Could not find taxonomy results for family {family}")
        if "errors" in node:
            for error in node["errors"]:
                logger.error(f"Error: {error['reason']}\n")
                sys.exit(100)

    return node["taxonomy"]["tax_id"]


def get_parent_taxid(taxid: int, ncbi_api_key: str | None):
    result = send_request_to_ncbi_taxonomy(taxid, ncbi_api_key)

    if "taxonomy_nodes" not in result:
        raise ValueError(f"Could not get taxonomy data for taxid {taxid}")
    node = result["taxonomy_nodes"][0]

    if "taxonomy" not in node:
        logger.info(f"Could not find taxonomy results for taxid {taxid}")
        if "errors" in node:
            for error in node["errors"]:
                logger.error(f"Error: {error['reason']}\n")
                sys.exit(100)

    if not node["taxonomy"].get("lineage"):
        logger.info(f"Could not find parent for taxid {taxid}")
        sys.exit(100)

    return node["taxonomy"]["lineage"][-1]


def get_assembly_stats(taxid: int, ncbi_api_key: str | None) -> dict:
    result = send_request_to_ncbi_genome_dataset_api(taxid, ncbi_api_key)

    assembly_stats = {}
    if result.get("reports") and isinstance(result["reports"], list):
        for report in result["reports"]:
            accession = report["accession"]
            report_assembly_stats = report.get("assembly_stats", {})
            assembly_stats[accession] = report_assembly_stats

    # if, for any reason, we could not get any assembly stat, we try with parent
    if not assembly_stats:
        logger.warning(
            f"Could not get assembly stats for children for taxid {taxid}. Trying with parent."
        )
        # recursive call with parent taxids
        parent_taxid = get_parent_taxid(taxid, ncbi_api_key)
        logger.info(f"Parent taxid: {parent_taxid}")
        return get_assembly_stats(parent_taxid, ncbi_api_key)

    # remove genbank assemblies whenever there is the equivalent refseq one
    filter_assembly_stats = {}
    for accession, report_assembl_stats in assembly_stats.items():
        if accession.startswith("GCA"):  # genbank
            refseq_accession = accession.replace("GCA", "GCF")
            if refseq_accession in assembly_stats:
                continue
        filter_assembly_stats[accession] = report_assembl_stats

    return filter_assembly_stats


def get_mean_assembly_length(assembly_stats: dict) -> float:
    assembly_lengths = [
        report_assembly_stats.get("total_sequence_length")
        for report_assembly_stats in assembly_stats.values()
        if report_assembly_stats.get("total_sequence_length") is not None
    ]
    assembly_lengths = [int(length) for length in assembly_lengths]
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

    family_taxid = get_family_taxid(family, args.ncbi_api_key)
    logger.info(f"Family taxid: {family_taxid}")

    logger.info(f"Getting assembly stats for children in taxid {family_taxid}")

    assembly_stats = get_assembly_stats(family_taxid, args.ncbi_api_key)

    mean_assembly_length = get_mean_assembly_length(assembly_stats)
    logger.info(f"Mean assembly length: {mean_assembly_length}")

    summary_dict = {
        "mean_assembly_length": mean_assembly_length,
    }

    output_dict = summary_dict | assembly_stats

    with open(args.outfile, "w") as fout:
        json.dump(output_dict, fout)

    logger.info("Done")
