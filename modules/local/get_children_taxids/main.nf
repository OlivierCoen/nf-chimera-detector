process GET_CHILDREN_TAXIDS {

    label 'process_single'

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/13675d7298d4faf2f44dd4dfb45106f4c2e3b49e98c0c167b3779f4a31195883/data':
        'community.wave.seqera.io/library/requests:2.32.4--9717c668df6a39b1' }"

    input:
    val family


    output:
    path("*.children_taxids.txt"),                                                                                    emit: taxids
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_all_children_taxids.py --family $family
    """

    stub:
    """
    touch test.children_taxids.txt
    """

}
