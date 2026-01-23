include { QUAST } from '../../../../modules/local/quast/main'

workflow RUN_QUAST {
    take:
    ch_main

    main:
    channel.empty().set { versions }
    /* prepare for quast:
     * This makes use of the input channel to obtain the reference and reference annotations
     * See quast module for details
     */
    channel.empty().set { quast_results }
    channel.empty().set { quast_tsv }

    ch_main
        .filter {
            it -> it.meta.quast
        }
        .multiMap { it ->
                quast_in: [
                    it.meta,
                    it.meta.qc_target,
                    it.meta.ref_fasta ?: [],
                    it.meta.ref_gff ?: [],
                    it.meta.ref_map_bam ?: [],
                    it.meta.assembly_map_bam
                ]
                use_ref: it.meta.use_ref
            }
        .set { quast_in }
    /*
    * Run QUAST
    */
    QUAST(quast_in.quast_in, quast_in.use_ref, false)
    QUAST.out.results.set { quast_results }
    QUAST.out.tsv.set { quast_tsv }
    QUAST.out.versions.set { versions }

    emit:
    quast_results
    quast_tsv
    versions
}
