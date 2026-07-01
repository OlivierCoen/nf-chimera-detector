process DOWNLOAD_ENA_FASTQ {

    label 'process_single'

    tag "${meta.family} :: txid${meta.taxid} :: ${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/9568ac2c728b11372f8d48c58ae686ac8802e3597145b161cbd50d9599605896/data':
        'community.wave.seqera.io/library/curl:8.21.0--ec591a19969f30d6' }"

    input:
    tuple val(meta), path(ena_ftp_url_file)

    output:
    tuple val(meta), path('*.fastq.gz'),                                                          emit: fastq
    tuple val("${task.process}"), val('curl'), eval("curl --version | head -1 | cut -d' ' -f2"), topic: versions

    script:
    """
    for url in \$(cat ${ena_ftp_url_file}); do
        echo "Downloading \${url}"
        wget \${url}
    done
    """

}
