include { RUN_MEDAKA } from '../run_medaka/main'
include { QC } from '../../../qc/main.nf'
include { RUN_LIFTOFF } from '../../../liftoff/main'

workflow POLISH_MEDAKA {
    take:
    ch_main
    meryl_kmers

    main:
    channel.empty().set { ch_versions }

    ch_main
        .filter {
            it -> it.meta.polish_medaka
        }
        .multiMap {
            it ->
            reads: [it.meta, it.meta.ontreads]
            reference: [it.meta, it.meta.assembly]
        }
        .set { ch_medaka_in }

    RUN_MEDAKA(ch_medaka_in.reads, ch_medaka_in.reference)

    RUN_MEDAKA.out.medaka_out.set { polished_assembly }

    polished_assembly
        .map { meta, polished_medaka -> [meta: meta + [ polished: [polished_medaka: polished_medaka ] ] ]}
        // After joining re-create the maps from the stored map
        .set { ch_medaka_out }

    ch_main
        .filter { it -> !it.polish_medaka }
        .mix(ch_medaka_out)
        .set { ch_main_out }

    ch_versions = ch_versions.mix(RUN_MEDAKA.out.versions)

    QC(
        ch_medaka_out.map { it -> [meta: it.meta - it.meta.subMap("assembly_map_bam") + [ assembly_map_bam: null] ] },
        polished_assembly.map { meta, polished -> [meta.id, polished] },
        meryl_kmers
    )


    ch_versions = ch_versions.mix(QC.out.versions)

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

    RUN_LIFTOFF(liftoff_in)

    ch_versions = ch_versions.mix(RUN_LIFTOFF.out.versions)

    versions = ch_versions

    emit:
    ch_main                 = ch_main_out
    quast_out               = QC.out.quast_out
    busco_out               = QC.out.busco_out
    merqury_report_files    = QC.out.merqury_report_files
    versions
}
