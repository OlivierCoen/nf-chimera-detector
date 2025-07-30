process SEQTK_SAMPLE {
    tag "${meta.taxid} :: ${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1aa6084536813ab32b6373b04407a898ab5de44692763b59ee0705e44974e0de/data' :
        'community.wave.seqera.io/library/seqtk_pigz:aa99a20f06d8e9a8' }"

    input:
    tuple val(meta), path(reads), val(sample_size)

    output:
    tuple val(meta), path("*.sampled.fastq.gz"),                                                            emit: reads
    tuple val("${task.process}"), val('seqtk'), eval("seqtk 2>&1 | awk 'NR==3' | sed 's/Version: //g'"),    topic: versions
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //g'"),              topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!(args ==~ /.*\ -s\ ?[0-9]+.*/)) {
        args += " -s100"
    }
    if ( !sample_size ) {
        error "SEQTK/SAMPLE must have a sample_size value included"
    }
    """
    printf "%s\\n" $reads | while read file;
    do
        FILE_STEM=\$(echo "\$file" | sed -E 's/\\.(fastq|fq)(\\.gz)?\$//')
        seqtk \\
            sample \\
            $args \\
            \$file \\
            $sample_size \\
            | pigz --no-name -p ${task.cpus} > \${FILE_STEM}.sampled.fastq.gz
    done
    """
}
