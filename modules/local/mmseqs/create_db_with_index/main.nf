process MMSEQS_CREATE_DB_WITH_INDEX {

    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1':
        'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1' }"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("${prefix}/"),                                                                           emit: db
    tuple val("${task.process}"), val('mmseqs'), eval("mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //'"),    topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = sequence.getExtension() == "gz" ? true : false
    def sequence_name = is_compressed ? sequence.getBaseName() : sequence
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${sequence} > ${sequence_name}
    fi

    mkdir -p ${prefix}

    mmseqs \\
        createdb \\
        ${args} \\
        ${sequence_name} \\
        ${prefix}/${prefix}

    mmseqs \\
        createindex \\
        ${args2} \\
        ${prefix}/${prefix} \\
        tmp \\
        --threads ${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    touch ${prefix}/${prefix}
    touch ${prefix}/${prefix}.dbtype
    touch ${prefix}/${prefix}.index
    touch ${prefix}/${prefix}.lookup
    touch ${prefix}/${prefix}.source
    touch ${prefix}/${prefix}_h
    touch ${prefix}/${prefix}_h.dbtype
    touch ${prefix}/${prefix}_h.index
    touch ${prefix}/${prefix}_h.idx
    """
}
