include { LIFTOFF } from '../../../modules/nf-core/liftoff/main'

workflow RUN_LIFTOFF {
  take:
  assembly
  inputs

  main:
  assembly
    .join(
      inputs.map { row -> [row.meta, row.ref_fasta, row.ref_gff] }
    )
    .set { liftoff_in }

  LIFTOFF(liftoff_in, [])

  LIFTOFF.out.gff3.set { lifted_annotations }

  emit:
  lifted_annotations
}
