process DOWNLOAD_NCBI_ASSEMBLY {

    label 'process_single'

    tag "${meta.taxid} :: $accession"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a6b13690259900baef6865722cb3a319103acc83b5bcab67504c88bde1e3a9f6/data':
        'community.wave.seqera.io/library/ncbi-datasets-cli_unzip:785aabe86637bae4' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path('*.{fasta,fa,fna}'),                                                                                   emit: assemblies
    tuple val("${task.process}"), val('ncbi-datasets-cli'), eval("datasets --version | sed 's/datasets version: //g'"),          topic: versions

    script:
    """
    datasets download genome accession $accession

    unzip -o ncbi_dataset.zip
    mv ncbi_dataset/data/${accession}/* .
    """

}
