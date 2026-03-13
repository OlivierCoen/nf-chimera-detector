#!/usr/bin/env bash

OUTFILE="failed_ids.txt"

sra_ids_string=$(nextflow log $(nextflow log | tail -1 | awk '{print $7}') -f name,status,workdir | grep "FAILED" | grep -oPe 'SR[RX]\d+')
sra_ids=($sra_ids_string)

srx_ids=()
for sra_id in "${sra_ids[@]}"
do
    if [[ "$sra_id" == "SRX"* ]]; then
        srx=$sra_id

    elif [[ "$sra_id" == "SRR"* ]]; then
        sra_url="https://www.ncbi.nlm.nih.gov/sra/?term=${sra_id}"
        tmpfile="${sra_id}.html"
        curl -s -o $tmpfile $sra_url
        srx=$(grep -oPe 'SRX\d+' $tmpfile | uniq)

    else
        echo "Unknown ID: ${sra_id}"
        exit 1
    fi

    srx_ids+=($srx)
done

printf "%s\n" "${srx_ids[@]}" > $OUTFILE
