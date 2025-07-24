process GET_BEST_NCBI_ASSEMBLY {

    label 'process_single'

    tag "$taxid"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8d/8da74f9621e013b83c1e2405c8e20fa68c625523b2336704784acaf81c724c04/data':
        'community.wave.seqera.io/library/jq_requests_tenacity:20a1d2f2027ac092' }"

    input:
    tuple val(meta), val(taxid)

    output:
    tuple val(meta), eval('jq -r .accession assembly_report.json'),                                                              emit: accession
    tuple val(meta), path('assembly_report.json'),                                                                               emit: assembly_report
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                                topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'),            topic: versions

    script:
    """
    get_best_assembly.py --taxon-id $taxid
    """

}
