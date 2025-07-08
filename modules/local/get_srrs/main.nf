process GET_SRRS {

    label 'process_single'

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/8770f64f21141a6537418e044fec40f635c9065828e726c93f2584c20f2ab09c/data':
        'community.wave.seqera.io/library/entrez-direct_pandas_tqdm:1af8e334fbe7f2f8' }"

    input:
    val taxid


    output:
    path("*_srrs.txt"),                                                                                               emit: srrs
    path("*_sra_metadata.csv"),                                                                                       emit: sra_metadata
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'),       topic: versions
    tuple val("${task.process}"), val('tqdm'),     eval('python3 -c "import tqdm; print(tqdm.__version__)"'),         topic: versions
    tuple val("${task.process}"), val('entrez-direct'), val('24.0'),                                                  topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_sra_metadata.py --taxon-id $taxid
    """

    stub:
    """
    touch test_srrs.txt test_sra_metadata.csv
    """

}
