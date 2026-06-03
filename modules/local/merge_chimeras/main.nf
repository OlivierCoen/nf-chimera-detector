process MERGE_CHIMERAS {

    label 'process_high'

    tag "${meta.family} :: txid${meta.taxid} :: ${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd1d206cd7587948d4761db523c452520f0e162468c9cf87055997e90c588fff/data':
        'community.wave.seqera.io/library/polars:1.41.1--128accb4274f01d4' }"

    input:
    tuple val(meta), path(files, stageAs: "**/*")

    output:
    tuple val(meta), path("${meta.id}.chimeras.csv"), emit: chimeras

    script:
    """
    merge_chimeras.py \\
        --chimeras ${files} \\
        --out ${meta.id}.chimeras.csv

    sleep 1
    """
}
