process GET_CHIMERA_READ_COVERAGE {

    label 'process_single'
    tag "${meta.id}"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b1b8abc03438fd8ccf1b79c2424f4f8fa73504e63d9d04b04dbb1dc200d8d2bf/data':
        'community.wave.seqera.io/library/matplotlib_pandas_polars_tqdm:f23204a31853358c' }"

    input:
    tuple val(meta), path(blast_hit_file), path(chimeras_file)

    output:
    tuple val(meta), path(prefix),                                                                                                             emit: coverage_dir
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                                               topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'),                                 topic: versions
    tuple val("${task.process}"), val('polars'), eval('python3 -c "import polars; print(polars.__version__)"'),                                 topic: versions
    tuple val("${task.process}"), val('matplotlib'), eval('python3 -c "import matplotlib; print(matplotlib.__version__)"'),                     topic: versions
    tuple val("${task.process}"), val('tqdm'), eval('python3 -c "import tqdm; print(tqdm.__version__)"'),                                       topic: versions

    script:
    def prefix = "${meta.id}"
    """
    get_chimera_read_coverage.py \\
        --hits $blast_hit_file \\
        --chimeras $chimeras_file \\
        --outdir $prefix
    """

    stub:
    """
    touch ${meta.id}/fake/test.raw.png
    touch ${meta.id}/fake/test.log.png
    """

}
