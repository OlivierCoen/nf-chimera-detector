process MERGE_CHIMERAS {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd1d206cd7587948d4761db523c452520f0e162468c9cf87055997e90c588fff/data':
        'community.wave.seqera.io/library/polars:1.41.1--128accb4274f01d4' }"

    input:
    tuple val(meta), path(files, stageAs: "**/*")

    output:
    tuple val(meta), path("${meta.id}.chimeras.csv"), emit: chimeras
    tuple val("${meta.family}"), val("${meta.id}"), env("NB_CHIMERAS"), topic: nb_chimeras

    script:
    """
    merge_chimeras.py \\
        --chimeras ${files} \\
        --out ${meta.id}.chimeras.csv

    NB_CHIMERAS=\$(cut -f 1 ${meta.id}.chimeras.csv | wc -l)
    """
}
