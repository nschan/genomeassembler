process LIFTOFF {
  tag "$meta.id"
  label 'process_high'
  
  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/liftoff:1.6.3--pyhdfd78af_0':
        'biocontainers/liftoff:1.6.3--pyhdfd78af_0' }"
  publishDir(
    path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
    mode: 'copy',
    overwrite: true,
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
  ) 
  conda "bioconda::liftoff=1.6.4"
  input:
      tuple val(meta), path(assembly), path(reference_fasta), path(reference_gff)
  
  output:
      tuple val(meta), path("*_liftoff.gff"), emit: lifted_annotations

  
  script:
      def prefix = task.ext.prefix ?: "${meta.id}"
  """
  if [[ ${assembly} == *.gz ]]; then
    zcat ${assembly} > assembly.fasta
  fi

  if [[ ${assembly} == *.fa || ${assembly} == *.fasta ]]; then
    cp ${assembly} assembly.fasta
  fi

  liftoff \\
    -g ${reference_gff} \\
    -p ${task.cpus} \\
    assembly.fasta  \\
    ${reference_fasta} \\
    -o ${assembly.baseName}_liftoff.gff
  """
}