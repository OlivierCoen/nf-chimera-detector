#!/usr/bin/env bash

set -u

SRA_ID=$1
ARGS=$2
MAX_ATTEMPTS=$3
DELAY=$4
MAX_TIME=$5

retry_with_backoff() {

    local attempt=1
    local output=
    local status=

    # Remove the first three arguments to this function in order to access
    # the 'real' command with `${@}`.
    shift 3

    while [ ${attempt} -le ${MAX_ATTEMPTS} ]; do
        output=$("${@}")
        status=${?}

        if [ ${status} -eq 0 ]; then
            break
        fi

        if [ ${attempt} -lt ${MAX_ATTEMPTS} ]; then
            echo "Failed attempt ${attempt} of ${MAX_ATTEMPTS}. Retrying in ${DELAY} s." >&2
            sleep ${DELAY}
        elif [ ${attempt} -eq ${MAX_ATTEMPTS} ]; then
            echo "Failed after ${attempt} attempts." >&2
            return ${status}
        fi

        attempt=$(( ${attempt} + 1 ))
        DELAY=$(( ${DELAY} * 2 ))
        if [ ${DELAY} -ge ${MAX_TIME} ]; then
            DELAY=${MAX_TIME}
        fi
    done

    echo "${output}"
}

export NCBI_SETTINGS="$PWD/!{ncbi_settings}"

retry_with_backoff prefetch $ARGS $SRA_ID

# check file integrity using vdb-validate or (when archive contains no checksums) md5sum
vdb-validate !{id} > vdb-validate_result.txt 2>&1 || exit 1
if grep -q "checksums missing" vdb-validate_result.txt; then
    VALID_MD5SUMS=$(curl --silent --fail --location --retry 3 --retry-delay 60 "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?filetype=run&acc=${SRA_ID}")
    LOCAL_MD5SUMS=$(md5sum ${SRA_ID}/* | cut -f1 -d' ')
    if ! grep -q -F -f <(echo "$LOCAL_MD5SUMS") <(echo "$VALID_MD5SUMS"); then
        echo "MD5 sum check failed" 1>&2
        exit 1
    fi
fi

