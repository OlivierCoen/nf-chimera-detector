process MEGAHIT {
    tag "${meta.taxid} :: ${meta.id}"
    label 'process_high'

    errorStrategy = {
        if (task.exitStatus == 100) {
            // ignoring cases when the assembly is empty
            log.warn("Assembly is empty for SRA ID ${meta.id}.")
            return 'ignore'
        }
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f2cb827988dca7067ff8096c37cb20bc841c878013da52ad47a50865d54efe83/data' :
        'community.wave.seqera.io/library/megahit_pigz:87a590163e594224' }"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("*.contigs.fa.gz"),                                                    emit: contigs
    tuple val("${task.process}"), val('megahit'), eval("megahit -v 2>&1 | sed 's/MEGAHIT v//'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end || !reads2 ? "-r ${reads1.join(',')}" : "-1 ${reads1.join(',')} -2 ${reads2.join(',')}"
    """
    megahit \\
        ${args} \\
        -t ${task.cpus} \\
        ${reads_command} \\
        --out-prefix ${prefix}

    # checking that the assembly is not empty
    [ -s megahit_out/${prefix}.contigs.fa ] || exit 100

    pigz \\
        --no-name \\
        -p ${task.cpus} \\
        ${args2} \\
        megahit_out/*.fa \\
        megahit_out/intermediate_contigs/*.fa

    mv megahit_out/* .
    """
}
