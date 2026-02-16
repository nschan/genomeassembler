process COUNT {
    tag "${meta.id}"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kmer-jellyfish:2.3.1--h4ac6f70_0'
        : 'biocontainers/kmer-jellyfish:2.3.1--h4ac6f70_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.jf"), emit: kmers
    tuple val("${task.process}"), val('jellyfish'), eval("jellyfish --version sed 's/jellyfish //'"), emit: versions_jellyfish, topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ ${fasta} == *.gz ]]; then
        zcat ${fasta} > ${fasta.baseName}.fasta
    fi
    if [[ ${fasta} == *.fa ]]; then
        cp ${fasta} ${fasta.baseName}.fasta
    fi
    if [[ ${fasta} == *.fastq ]]; then
        cp ${fasta} ${fasta.baseName}.fasta
    fi
    jellyfish count \\
        -m ${meta.jellyfish_k} \\
        -s 140M \\
        -C \\
        -t ${task.cpus} ${fasta.baseName}.fasta
    mv mer_counts.jf ${prefix}_mer_counts.jf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mer_counts.jf
    """

}
