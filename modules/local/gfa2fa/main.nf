process GFA_2_FA {
    tag "${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(gfa_file)

    output:
    tuple val(meta), path("*fa.gz"), emit: contigs_fasta
    tuple val("${task.process}"), val('awk'), eval("awk |& head -n1 | grep -o 'v[0-9\\.]*' | sed 's/v//'"), emit: versions_awk, topic: versions
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | head -n1 | sed 's/bgzip (htslib) //'"), emit: versions_gzip, topic: versions

    script:
    """
    outfile=\$(basename $gfa_file .gfa).fa.gz
    awk '/^S/{print ">"\$2;print \$3}' ${gfa_file} \\
    | bgzip > \$outfile
    """

    stub:
    """
    outfile=\$(basename $gfa_file .gfa).fa.gz
    touch \$outfile
    """
}
