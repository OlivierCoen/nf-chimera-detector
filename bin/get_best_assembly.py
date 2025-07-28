#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import requests
import argparse
import json
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
NCBI_GENOME_DATASET_REPORT_API_PARAMS = "page_size=1000"

NCBI_API_HEADERS = {
    "accept": "application/json",
    "content-type": "application/json"
}


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get best assembly for a specific taxon ID"
    )
    parser.add_argument(
        "--taxon-id", type=int, dest="taxid", required=True, help="Taxon ID on NCBI Taxonomy"
    )
    parser.add_argument(
        "--out", dest="outfile", type=str, required=True, help="Outfile name"
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


def get_assembly_with_best_stats(reports: list[dict]):
    sorted_reports = sorted(
        reports,
        key=lambda x: (
            int(x.get('assembly_stats').get('total_sequence_length', 0)),
            -int(x.get('assembly_stats', {}).get('total_number_of_chromosomes', 1E9)),
        ),
        reverse=True
    )
    return sorted_reports[0]


def get_current_assemblies(reports: list[dict]):
    current_assembly_reports = [
        report for report in reports
        if report.get('assembly_info', {}).get('assembly_status') == "current"
    ]
    if not current_assembly_reports:
        return reports

    refseq_reports = [
        report for report in current_assembly_reports
        if report.get('source_database') == "SOURCE_DATABASE_REFSEQ"
    ]

    if refseq_reports:
        return refseq_reports
    else:
        return current_assembly_reports


def get_reference_assembly(reports: list[dict]):
    current_assembly_report = get_current_assemblies(reports)
    return get_assembly_with_best_stats(current_assembly_report)

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()
    taxid = args.taxid

    logger.info(f"Getting best NCBI assembly for taxid: {taxid}")
    result = send_request_to_ncbi_genome_dataset_api(taxid)

    try:
        reports = result['reports']
        reference_assembly_report = get_reference_assembly(reports)
        logger.info(f"Best assembly: {reference_assembly_report['accession']}")
    except Exception as e:
        logger.error(f"Could not get any assembly for taxid {taxid}: {e}")
        reference_assembly_report = {"accession": "NONE"}

    with open(args.outfile, 'w') as fout:
        json.dump(reference_assembly_report, fout)

    logger.info("Done")
