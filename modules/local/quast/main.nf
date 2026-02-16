process QUAST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a515d04307ea3e0178af75132105cd36c87d0116c6f9daecf81650b973e870fd/data' :
        'community.wave.seqera.io/library/quast:5.3.0--755a216045b6dbdd' }"

    input:
    tuple val(meta), path(consensus), path(fasta), path(gff), path(ref_bam), path(bam)
    val use_fasta
    val use_gff

    output:
    path "${meta.id}*/*", emit: results
    path "*report.tsv", emit: tsv
    tuple val("${task.process}"), val('quast'), eval("quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//' | tail -n1"), emit: versions_medaka, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def features = use_gff ? "--features ${gff}" : ''
    def reference = use_fasta ? "-r ${fasta}" : ''
    def reference_bam = ref_bam ? "--ref-bam ${ref_bam}" : ''

    """
    quast.py \\
        --output-dir ${prefix} \\
        ${reference} \\
        ${features} \\
        --threads ${task.cpus} \\
        ${consensus.join(' ')} \\
        --glimmer \\
        ${reference_bam} \\
        --bam ${bam} \\
        --large \\
        ${args}

    ln -s ${prefix}/report.tsv ${prefix}_report.tsv
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix} && touch ${prefix}/report.tsv
    ln -s ${prefix}/report.tsv ${prefix}_report.tsv
    """
}
