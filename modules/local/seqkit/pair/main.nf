process SEQKIT_PAIR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/85b40b925e4d4a62f9b833bbb0646d7ea6cf53d8a875e3055f90da757d7ccd27/data':
        'community.wave.seqera.io/library/seqkit:2.12.0--5e60eb93d3a59212' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("result/*"),                                                      emit: reads
    tuple val("${task.process}"), val('seqkit'), eval("seqkit | sed '3!d; s/Version: //'"), topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    seqkit \\
        pair \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --out-dir result \\
        --threads ${task.cpus} \\
        ${args}

    """
}
