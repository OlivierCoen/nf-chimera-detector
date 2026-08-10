process DOWNLOAD_ENA_FASTQ {

    label 'process_high'
    tag "${meta.family} :: txid${meta.taxid} :: ${meta.sra_id}"

    maxForks 1

    errorStrategy {
        if ( task.exitStatus == 1 ) {
            'retry'
        } else if ( task.exitStatus in ( [104, 175] + (130..145).toList() ) ) { // OOM & related errors; should be retried as long as memory does not fit
            'retry'
        } else {
            'ignore'
        }
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e19ee5cc46963e10871d03bbcc4d41823f4a57ca6f52d7bcdf3147d0e895d8a/data':
        'community.wave.seqera.io/library/axel:2.17.13--3a206f517a443ff3' }"

    input:
    tuple val(meta), path(ena_ftp_url_file)

    output:
    tuple val(meta), path('*.fastq.gz'),                                                          emit: fastq
    tuple val("${task.process}"), val('axel'), eval("axel --version | head -1 | cut -d' ' -f2"), topic: versions

    script:
    """
    for url in \$(cat ${ena_ftp_url_file}); do
        # concert ftp URL to https to avoid ftp connection issues
        #http_url=\$(echo \$url | sed 's#ftp://#https://#g')
        echo "Downloading \${http_url}"
        axel \\
            -n ${task.cpus} \\
            \${http_url}
    done
    """

}
