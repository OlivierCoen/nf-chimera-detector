process NCBI_ASSEMBLY_STATS {

    label 'process_single'

    tag "$family"

    errorStrategy = {
        if (task.exitStatus == 100) {
            // ignoring cases when family does not have children
            log.warn("Could not find species for family ${family}.")
            return 'ignore'
        }
    }

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8d/8da74f9621e013b83c1e2405c8e20fa68c625523b2336704784acaf81c724c04/data':
        'community.wave.seqera.io/library/jq_requests_tenacity:20a1d2f2027ac092' }"

    input:
    val family

    output:
    tuple val(family), eval('jq -r .mean_assembly_length *.stats.json'),                                                                        emit: mean_lengths
    tuple val(family), path('*.stats.json'),                                                                                                    emit: stats
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                                               topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'),                           topic: versions
    tuple val("${task.process}"), val('tenacity'), eval('python3 -c "from importlib.metadata import version; print(version(\'tenacity\'))"'),   topic: versions
    tuple val("${task.process}"), val('jq'),       eval("jq --version | sed 's/jq-//g'"),                                                       topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "$family"
    """
    get_assembly_stats.py --family $family --out ${prefix}.stats.json
    """
}
