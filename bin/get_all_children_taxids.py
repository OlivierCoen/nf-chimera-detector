#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import requests
import subprocess
import sys
import argparse
import logging

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

NCBI_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
NCBI_API_HEADERS = {
    "accept": "application/json",
    "content-type": "application/json"
}
CHUNKSIZE = 2000


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
        "--family", type=str, required=True, help="Family name"
    )
    return parser.parse_args()


def get_virus_taxids(taxid: str):
    command = f'esearch -db taxonomy -query "txid{taxid}[Subtree]" | efetch -format uid'
    result = subprocess.run(command, shell=True, capture_output=True, check=True, text=True)
    if result.returncode != 0:
        logger.info(result.stderr)
        raise RuntimeError('Could not get virus taxids')
    return [taxid for taxid in result.stdout.split('\n') if taxid != '']


def send_request_to_ncbi_taxonomy(taxons: list[str]):
    data = {
        "taxons": taxons
    }
    response = requests.post(NCBI_API_URL, headers=NCBI_API_HEADERS, json=data)
    # check status
    if not response.ok:
        raise RuntimeError(f"Error: {response.status_code}. {response.text}")
    return response.json()


def get_children_taxids(taxids: list[str]):
    # querying NCBI taxonomy API
    result = send_request_to_ncbi_taxonomy(taxids)

    # parsing result
    children_taxids = []
    for taxonomy_node_dict in result.get("taxonomy_nodes", []):
        if (taxon_dict := taxonomy_node_dict.get("taxonomy")) is not None:

            taxid = taxon_dict.get("tax_id")
            rank = taxon_dict.get("rank")

            if taxid is not None and rank == 'SPECIES':
                children_taxids.append(taxid)

    return children_taxids


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()
    
    family = args.family

    logger.info(f"Getting all children taxids for family {family}")
    family_dna_children_taxids = []

    result = send_request_to_ncbi_taxonomy([family])
    if len(result['taxonomy_nodes']) > 1:
        raise ValueError(f"Multiple taxids for family {family}")
    node = result['taxonomy_nodes'][0]

    if "taxonomy" not in node:
        logger.info(f"Could not find taxnomomy results for family {family}")
        if "errors" in node:
            for error in node["errors"]:
                logger.error(f"Error: {error['reason']}\n")
                sys.exit(100)

    family_taxid = node['taxonomy']['tax_id']
    logger.info(f"Family taxid: {family_taxid}")
    virus_taxids = get_virus_taxids(family_taxid)

    virus_taxids_chunks = [virus_taxids[i:i + CHUNKSIZE] for i in range(0, len(virus_taxids), CHUNKSIZE)]
    for virus_taxids_chunk in virus_taxids_chunks:
        family_dna_children_taxids += get_children_taxids(virus_taxids_chunk)

    logger.info(f"Obtained {len(family_dna_children_taxids)} children taxids\n")

    outfile = f"{family}.children_taxids.txt"
    with open(outfile, 'w') as fout:
        for taxid in family_dna_children_taxids:
            fout.write(f"{taxid}\n")

    logger.info("Done")
