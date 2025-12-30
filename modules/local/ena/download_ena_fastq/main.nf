process DOWNLOAD_ENA_FASTQ {

    label 'process_single'

    tag "${meta.family} :: txid${meta.taxid} :: ${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4cac47b5cb5de60c6e0f069ee809ac1c08b230fac691e145c106c3faa747129d/data':
        'community.wave.seqera.io/library/aria2:1.34.0--555114d61d424337' }"

    input:
    tuple val(meta), path(ena_ftp_url_file)

    output:
    tuple val(meta), path('*.fastq.gz'),                                                              emit: fastq
    tuple val("${task.process}"), val('aria2'), eval('aria2c --version | head -1 | cut -d" " -f3'),   topic: versions

    script:
    """
    for url in \$(cat ${ena_ftp_url_file}); do
        echo "Downloading \${url}"
        aria2c \\
            -x ${task.cpus} \\
            -s ${task.cpus} \\
            -o \${url##*/} \\
            \${url}
    done
    """

}
