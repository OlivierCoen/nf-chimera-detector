#!/usr/bin/env bash

VIRUSES_TAXON_ID=$(datasets summary taxonomy taxon viruses | jq '.reports[0].taxonomy.classification.acellular_root.id')

