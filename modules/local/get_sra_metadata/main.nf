process GET_SRA_METADATA {

    label 'process_single'

    tag "$taxid"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/7140956576e6779ff2bea610b8e41bde95d4909d67b31634d0ddb6dba50aef5a/data':
        'community.wave.seqera.io/library/requests_tenacity_xmltodict:9e74a2aeeb88aab9' }"

    input:
    val taxid


    output:
    path("*.sra_ids.txt"),                                                                                                       emit: sra_id_files
    path("*.sra_metadata.json"),                                                                                                 emit: sra_metadata
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                                topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'),            topic: versions
    tuple val("${task.process}"), val('xmltodict'), eval('python3 -c "import xmltodict; print(xmltodict.__version__)"'),         topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_sra_metadata.py --taxon-id $taxid
    """

    stub:
    """
    touch test.sra_ids.txt test.sra_metadata.csv
    """

}
