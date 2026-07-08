process SEQKIT_PAIR {
    tag "$meta.id"
    label 'process_medium'

    errorStrategy {
        if ( task.exitStatus == 255 && task.attempt <= 2 ) {
            'retry'
        } else if ( task.exitStatus in ( [104, 175] + (130..145).toList() ) ) { // OOM & related errors; should be retried as long as memory does not fit
            'retry'
        } else {
            'ignore'
        }
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/85b40b925e4d4a62f9b833bbb0646d7ea6cf53d8a875e3055f90da757d7ccd27/data':
        'community.wave.seqera.io/library/seqkit:2.12.0--5e60eb93d3a59212' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("result/*"),                                                      emit: reads
    tuple val("${task.process}"), val('seqkit'), eval("seqkit | sed '3!d; s/Version: //'"), topic: versions

    script:
    def args = task.ext.args ?: ''
    def sanitise_first = task.exitStatus == 255 ? 'true' : 'false'
    if ( sanitise_first == 'true' ) {
        // getting full file extension
        // if file extension is gz, get like 'fastq.gz', else get like 'fastq'
        def complete_suffix = reads[0].getExtension() == "gz" ? reads[0].baseName.tokenize('.')[-1] + ".gz" : reads[0].getExtension()
        """
        echo "Sanitising ${reads[0]} and ${reads[1]}"
        mkdir -p sanitised
        for read_file in ${reads[0]} ${reads[1]}
        do
            outfile="sanitised/\$(basename \$read_file)"
            seqkit \\
                sana \\
                \$read_file \\
                --out-file \$outfile \\
                --threads ${task.cpus} \\
                ${args}
        done

        seqkit \\
            pair \\
            -1 sanitised/${reads[0]} \\
            -2 sanitised/${reads[1]} \\
            --out-dir result \\
            --threads ${task.cpus} \\
            ${args}

        echo "Removing sanitised files"
        rm -rf sanitised
        """
    } else {
        """
        seqkit \\
            pair \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --out-dir result \\
            --threads ${task.cpus} \\
            ${args}
        """
    }
}
