process FIND_CHIMERAS {

    label 'process_low'

    tag "${meta.taxid} :: ${meta.id}"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/10eb3d3e0d620ea2ca74fa1bbb415026191ace60e049b7f264293acc2bd93818/data':
        'community.wave.seqera.io/library/r-base_r-data.table_r-optparse:574927877459455e' }"

    input:
    tuple val(meta), path("blast_hits.against_target.txt"), path("blast_hits.against_assembly.txt")

    output:
    tuple val(meta), path("*_chimeras.csv"),                                                                                  emit: csv, optional: true
    tuple val("${task.process}"), val('R'),          eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),       topic: versions
    tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_chimeras"
    """
    # checking that both files are not empty

    find_chimeras.R \\
        --target-hits blast_hits.against_target.txt \\
        --assembly-hits blast_hits.against_assembly.txt \\
        --out ${prefix}.csv
    """



}
