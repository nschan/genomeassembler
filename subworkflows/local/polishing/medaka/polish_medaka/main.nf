include { MEDAKA_PARALLEL as MEDAKA } from '../../../../../modules/local/medaka/medaka_consensus/main'
include { QC } from '../../../qc/main.nf'
include { LIFTOFF } from '../../../../../modules/nf-core/liftoff/main'

workflow POLISH_MEDAKA {
    take:
    ch_main
    meryl_kmers

    main:
    channel.empty().set { ch_versions }

    ch_main
        .map {
            it ->
            [ it.meta, it.meta.ontreads, it.meta.assembly ]

        }
        .set { ch_medaka_in }

    MEDAKA(ch_medaka_in)

    MEDAKA.out.assembly.set { polished_assembly }

    polished_assembly
        .map { meta, polished_medaka -> [meta: meta + [ polished: [medaka: polished_medaka ] ] ]}
        // After joining re-create the maps from the stored map
        .set { ch_medaka_out }

    ch_medaka_out
        .set { ch_main_out }

    QC(
        ch_medaka_out.map { it -> [meta: it.meta - it.meta.subMap("assembly_map_bam") + [ assembly_map_bam: null] ] },
        polished_assembly.map { meta, polished -> [meta.id, polished] },
        meryl_kmers
    )


    ch_medaka_out
        .filter {
            it -> it.meta.lift_annotations
        }
        .map { it ->
                [
                it.meta,
                it.meta.polished.medaka,
                it.meta.ref_fasta,
                it.meta.ref_gff
                ]
        }
        .set { liftoff_in }

    LIFTOFF(liftoff_in, [])

    emit:
    ch_main                 = ch_main_out
    quast_out               = QC.out.quast_out
    busco_out               = QC.out.busco_out
    merqury_report_files    = QC.out.merqury_report_files
}
