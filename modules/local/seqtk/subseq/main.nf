process SEQTK_SUBSEQ {
    tag "$sequences"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1aa6084536813ab32b6373b04407a898ab5de44692763b59ee0705e44974e0de/data' :
        'community.wave.seqera.io/library/seqtk_pigz:aa99a20f06d8e9a8' }"

    input:
    tuple val(meta), path(sequences), path(filter_list)

    output:
    tuple val(meta), path("*.gz"),                                                                          emit: sequences
    tuple val("${task.process}"), val('seqtk'), eval("seqtk 2>&1 | awk 'NR==3 {print \$2}'"),               topic: versions
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //g'"),              topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${sequences.simpleName}_filtered"
    def extension = file( sequences.baseName ).getExtension()
    """
    seqtk \\
        subseq \\
        $args \\
        $sequences \\
        $filter_list | \\
        pigz --no-name > ${prefix}.${extension}.gz
    """
}
