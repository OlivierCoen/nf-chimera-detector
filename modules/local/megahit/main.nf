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
    tuple val("${meta.family}"), val("${meta.id}"), env('GENOME_NB_BASES'),                      topic: assembled_genome_nb_bases
    path("*.log"),                                                                               topic: megahit_multiqc
    tuple val("${task.process}"), val('megahit'), eval("megahit -v 2>&1 | sed 's/MEGAHIT v//'"), topic: versions

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

    assembly=megahit_out/${prefix}.contigs.fa

    # checking that the assembly is not empty
    [ -s \$assembly ] || exit 100

    # compute total genome size
    GENOME_NB_BASES=\$(grep -v "^>" \$assembly | tr -d "\n" | wc -m)

    pigz \\
        --no-name \\
        -p ${task.cpus} \\
        ${args2} \\
        megahit_out/*.fa \\
        megahit_out/intermediate_contigs/*.fa

    mv megahit_out/* .
    """
}
