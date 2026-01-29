process DORADO_ALIGNER {
    tag "${meta.id}"
    label 'process_high'

    container "docker.io/nanoporetech/dorado:shaf2aed69855de85e60b363c9be39558ef469ec365"

    input:
    tuple val(meta), path(ref), path(reads)

    output:
    tuple val(meta), path("${meta.id}_dorado_aligned.bam"), emit: bam
    tuple val(meta), path("${meta.id}_dorado_aligned.bai"), emit: bai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    dorado aligner \\
        -t ${task.cpus} \\
        ${ref} \\
        ${reads} \\
        ${args} \\
        > ${meta.id}_dorado_aligned.bam
    samtools index ${meta.id}_dorado_aligned.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: "\$(dorado --version 2>&1 | head -n1)"
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}/${prefix}.bam
    touch ${prefix}/${prefix}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: "\$(dorado --version 2>&1 | head -n1)"
    END_VERSIONS
    """
}
