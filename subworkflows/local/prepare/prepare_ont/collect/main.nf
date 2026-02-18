include { COLLECT_READS } from '../../../../../modules/local/collect_reads/main'

workflow COLLECT {
    take:
    ch_input

    main:
    ch_input
        .filter {
            it -> it.ont_collect
        }
        .map { row -> [row.meta, row.meta.ontreads] }
        .set { reads }

    COLLECT_READS(reads)
    COLLECT_READS.out.combined_reads.set { reads }

    emit:
    reads
}
