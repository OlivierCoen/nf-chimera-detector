process PARSE_MMSEQS_OUTPUT {

    label 'process_single'

    tag "${meta.family} :: ${meta.id}"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c03d8fbb44376692e10bb9e56be2577fc35446c193637735de9fed182e6b58df/data':
        'community.wave.seqera.io/library/pandas:2.3.1--139e2fa6c1f18206' }"

    input:
    tuple val(meta), path(file)
    val max_evalue

    output:
    tuple val(meta), eval('cat status.txt'),                                                                           emit: status
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                          topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'),            topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_mmseqs_output.py \\
        --file $file \\
        --max-evalue $max_evalue
    """

    stub:
    """
    echo "PASS" >> status.txt
    """

}
