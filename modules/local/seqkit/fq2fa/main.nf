process SEQKIT_FQ2FA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/85b40b925e4d4a62f9b833bbb0646d7ea6cf53d8a875e3055f90da757d7ccd27/data' :
        'community.wave.seqera.io/library/seqkit:2.12.0--5e60eb93d3a59212' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), env("READ_FASTA_SUM_LEN"), path("*.fa.gz"),                            emit: fasta
    tuple val("${task.process}"), val('seqkit'), eval("seqkit | sed '3!d; s/Version: //'"), topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqkit \\
        fq2fa \\
        $args \\
        -j $task.cpus \\
        -o ${prefix}.fa.gz \\
        $fastq

    READ_FASTA_SUM_LEN=\$(seqkit stats *.fa.gz | tail -1 | tr -s '[:space:]' '\t' | cut -f5 | sed 's/,//g')
    """

    stub:
    def prefix = task.ext.prefix ?: "${fastq.simpleName}"
    """
    echo "" | gzip > ${prefix}.fa.gz
    """
}
