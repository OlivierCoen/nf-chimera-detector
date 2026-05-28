process MERGE_CHIMERAS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(files, stageAs: "**/*")

    output:
    tuple val(meta), path("${meta.id}.chimeras.csv"), emit: chimeras
    tuple val("${meta.family}"), val("${meta.id}"), env("NB_CHIMERAS"), topic: nb_chimeras

    script:
    """
    cat ${files} > ${meta.id}.chimeras.csv
    NB_CHIMERAS=\$(cut -f 1 ${meta.id}.chimeras.csv | wc -l)
    """
}
