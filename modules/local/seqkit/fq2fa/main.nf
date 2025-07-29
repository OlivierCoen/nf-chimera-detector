process SEQKIT_FQ2FA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fa.gz"),                                                     emit: fasta
    tuple val("${task.process}"), val('pigz'), eval("seqkit | sed '3!d; s/Version: //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fastq.simpleName}"

    """
    seqkit \\
        fq2fa \\
        $args \\
        -j $task.cpus \\
        -o ${prefix}.fa.gz \\
        $fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${fastq.simpleName}"
    """
    echo "" | gzip > ${prefix}.fa.gz
    """
}
