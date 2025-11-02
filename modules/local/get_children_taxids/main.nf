process GET_CHILDREN_TAXIDS {

    label 'process_single'

    tag "$family"

    errorStrategy {
        if (task.exitStatus == 100) {
            // ignoring cases when family does not have children
            log.warn("Could not find species for family ${family}.")
            return 'ignore'
        }
    }

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b4d686ef63e22bc4d461178fc241cefddd2aa3436e189d3787c8e019448f056e/data':
        'community.wave.seqera.io/library/requests_tenacity_tqdm:126dbed8ef3ff96f' }"

    input:
    tuple val(meta), val(family)

    output:
    tuple val(meta), path("*.taxids2names.csv"), emit: taxid_to_names_files

    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                                               topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'),                           topic: versions
    tuple val("${task.process}"), val('tenacity'), eval('python3 -c "from importlib.metadata import version; print(version(\'tenacity\'))"'),   topic: versions

    script:
    """
    get_children_taxids.py --family $family
    """

    stub:
    """
    touch test.species_taxids.txt
    touch test.taxids2names.csv
    """

}
