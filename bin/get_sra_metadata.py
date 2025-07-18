#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import json

import requests
import sys
import xml.etree.ElementTree as ET
import xmltodict
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

# E-UTILITIES OL API
ESEARCH_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={db}&id={id}"
ESEARCH_RETMAX = 1000000000 # max retmax that worked
CHUNKSIZE = 2000

SRA_IDS_OUTFILE_SUFFIX = ".sra_ids.txt"
EXPERIMENT_OUTFILE_SUFFIX = ".sra_metadata.json"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute general statistics from count data for each sample"
    )
    parser.add_argument(
        "--taxon-id", type=int, dest="taxid", required=True, help="Taxon ID on NCBI Taxonomy"
    )
    return parser.parse_args()

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
    params = dict(
        db=database,
        term=query,
        retmax=ESEARCH_RETMAX
    )
    response = requests.get(ESEARCH_BASE_URL, params=params)
    response.raise_for_status()
    return response.text


@retry(
    retry=retry_if_exception_type(requests.exceptions.HTTPError),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_efetch_query(query_ids: list[str], database: str):
    """
    Query NCBI's db with Efetch API
    """
    query_ids_str = ",".join(query_ids)
    url = EFETCH_BASE_URL.format(db=database, id=query_ids_str)
    response = requests.get(url)
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
    return [id_element.text for id_element in root.findall('.//Id')]


def parse_sra_accessions_from_xml(xml_string: str) -> list[dict]:
    xml_dict = xmltodict.parse(xml_string)
    experiments = xml_dict['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
    return [
        experiment_dict['EXPERIMENT'] for experiment_dict in experiments
        if isinstance(experiment_dict, dict)
    ]



def get_sra_ids_from_taxid(taxid: str, dna_only: bool) -> list[str]:
    """
    Get list of SRA experiment IDs given a NCBI taxonomy ID
    :param taxid:
    :return: list of SRA experiment IDs
    """
    query = f'txid{taxid}'
    if dna_only:
        query += ' AND "biomol dna"[Properties]'
    xml_string = send_esearch_query(
        query=query,
        database="sra"
    )
    return parse_ids_from_xml(xml_string)


def fetch_sra_experiments(sra_uids: list[str]) -> list[dict]:
    """
    Get list of SRA experiment IDs given a NCBI taxonomy ID
    :param taxid:
    :return: list of SRA experiment IDs
    """
    xml_string = send_efetch_query(
        query_ids=sra_uids,
        database="sra"
    )
    return parse_sra_accessions_from_xml(xml_string)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    logger.info(f"Fetching sra experiment UIDs from taxon ID {args.taxid}")
    sra_uids = get_sra_ids_from_taxid(args.taxid, dna_only=True)
    logger.info(f"Got {len(sra_uids)} SRA experiment UIDs")

    logger.info("Fetching sra experiment metadata for each SRA experiment ID")
    experiments = []
    sra_uids_chunks = [sra_uids[i:i + CHUNKSIZE] for i in range(0, len(sra_uids), CHUNKSIZE)]
    for sra_uids_chunk in sra_uids_chunks:
        experiments += fetch_sra_experiments(sra_uids_chunk)

    outfile = f"{args.taxid}{EXPERIMENT_OUTFILE_SUFFIX}"
    with open(outfile, 'w') as fout:
        json.dump(experiments, fout)

    srrs = [exp['@accession'] for exp in experiments]

    outfile = f"{args.taxid}{SRA_IDS_OUTFILE_SUFFIX}"
    with open(outfile, 'w') as fout:
        for srr in srrs:
            fout.write(f"{srr}\n")

    logger.info("Done")
