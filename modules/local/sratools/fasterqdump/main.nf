process SRATOOLS_FASTERQDUMP {
    tag "${meta.taxid} :: ${sra.name}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37aacd127aa32161d8b38a83efb18df01a8ab1d769a93e88f80342d27801b548/data' :
        'community.wave.seqera.io/library/sra-tools_pigz:4a694d823f6f7fcf' }"

    input:
    tuple val(meta), path(sra)
    path ncbi_settings

    output:
    tuple val(meta), path('*.fastq.gz'),                                                                        emit: reads
    tuple val("${task.process}"), val('sratools'), eval("fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+'"),    topic: versions
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //g'"),                  topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${sra.name}"
    """
    export NCBI_SETTINGS="\$PWD/${ncbi_settings}"

    fasterq-dump \\
        $args \\
        --split-files \\
        --threads $task.cpus \\
        --outfile ${prefix}.fastq \\
        ${sra}

    pigz \\
        $args2 \\
        --no-name \\
        --processes $task.cpus \\
        *.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${sra.name}"
    """
    touch ${prefix}.fastq.gz
    """
}
