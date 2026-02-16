include { LIFTOFF } from '../../../modules/nf-core/liftoff/main'

workflow RUN_LIFTOFF {
    take:
    liftoff_in

    main:
    LIFTOFF(liftoff_in, [])

    LIFTOFF.out.gff3.set { lifted_annotations }

    emit:
    lifted_annotations
}
