process GET_CHIMERA_READ_COVERAGE {

    label 'process_high_memory'
    tag "${meta.id}"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/82/8285150b31c413e715428d9b65b13b277941be75f4923bfc4388d202ad0a7bcf/data':
        'community.wave.seqera.io/library/matplotlib_pandas_polars_pyarrow_tqdm:a426e1ee01722428' }"

    input:
    tuple val(meta), path(blast_hit_file), path(chimeras_file)

    output:
    tuple val(meta), path("*"), optional: true,                                                                              emit: coverage_dir
    tuple val("${task.process}"), val('python'),     eval("python3 --version | sed 's/Python //'"),                          topic: versions
    tuple val("${task.process}"), val('pandas'),     eval('python3 -c "import pandas; print(pandas.__version__)"'),          topic: versions
    tuple val("${task.process}"), val('polars'),     eval('python3 -c "import polars; print(polars.__version__)"'),          topic: versions
    tuple val("${task.process}"), val('matplotlib'), eval('python3 -c "import matplotlib; print(matplotlib.__version__)"'),  topic: versions
    tuple val("${task.process}"), val('tqdm'),       eval('python3 -c "import tqdm; print(tqdm.__version__)"'),              topic: versions

    script:
    def prefix = ""
    """
    export MPLCONFIGDIR=\${PWD}

    get_chimera_read_coverage.py \\
        --hits $blast_hit_file \\
        --chimeras $chimeras_file \\
        --srr ${meta.id}
    """

    stub:
    """
    touch ${meta.id}/fake/test.raw.png
    touch ${meta.id}/fake/test.log.png
    """

}
