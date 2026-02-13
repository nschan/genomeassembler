include { BUSCO_BUSCO as BUSCO } from '../../../../modules/nf-core/busco/busco/main'

workflow RUN_BUSCO {
    take:
    ch_main

    main:
    channel.empty().set { batch_summary }
    channel.empty().set { short_summary_txt }
    channel.empty().set { short_summary_json }

    ch_main
        .filter {
            it -> it.meta.busco
        }
        .multiMap { it ->
                fasta: [
                    it.meta,
                    it.meta.qc_target
                ]
                busco_lineage: it.meta.busco_lineage
                busco_db: it.meta.busco_db ? file(it.meta.busco_db, checkIfExists: true) : []
            }
        .set { busco_in }

    BUSCO(busco_in.fasta, 'genome', busco_in.busco_lineage, busco_in.busco_db , [], true)
    BUSCO.out.batch_summary.set { batch_summary }
    BUSCO.out.short_summaries_txt.set { short_summary_txt }
    BUSCO.out.short_summaries_json.set { short_summary_json }

    emit:
    batch_summary
    short_summary_json
    short_summary_txt
}
