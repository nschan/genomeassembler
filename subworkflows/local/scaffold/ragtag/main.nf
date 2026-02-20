include { RAGTAG_SCAFFOLD } from '../../../../modules/nf-core/ragtag/scaffold/main'
include { QC } from '../../qc/main'
include { LIFTOFF } from '../../../../modules/nf-core/liftoff/main'


workflow RUN_RAGTAG {
    take:
    ch_main
    meryl_kmers

    main:
    ch_main
        .multiMap { it ->
                    def assembly_to_scaffold =
                                it.meta.scaffold ?
                                (
                                    it.meta.scaffolds_hic ?:
                                    it.meta.scaffolds_longstitch ?:
                                    it.meta.scaffolds_links
                                ) :
                                it.meta.polished ?
                                (
                                    it.meta.polished.pilon ?:
                                    it.meta.polished.medaka ?:
                                    it.meta.polished.dorado
                                ) :
                                it.meta.assembly
                    assembly:
                        [
                            it.meta,
                            assembly_to_scaffold
                        ]
                    reference:
                        [
                            it.meta,
                            it.meta.ref_fasta
                        ]
                    }
        .set { ragtag_in }

    ragtag_in.assembly.dump(tag: "SCAFFOLD: RAGTAG: Assembly inputs")
    ragtag_in.reference.dump(tag: "SCAFFOLD: RAGTAG: Reference inputs")

    RAGTAG_SCAFFOLD(ragtag_in.assembly, ragtag_in.reference, [[], []], [[], [], []])

    RAGTAG_SCAFFOLD.out.corrected_assembly
        .map { meta, corrected -> [meta: meta + [ scaffolds_ragtag: corrected] ] }
        .set { ch_main_scaffolded }

    QC(ch_main_scaffolded.map { it -> [meta: it.meta - it.meta.subMap("assembly_map_bam") + [assembly_map_bam: null] ] },
        RAGTAG_SCAFFOLD.out.corrected_assembly.map { meta, corrected -> [ meta.id, corrected ] },
        meryl_kmers)


    ch_main_scaffolded
        .filter {
            it -> it.lift_annotations
        }
        .map { it ->
                [
                it.meta,
                it.meta.scaffolds_ragtag,
                it.meta.ref_fasta,
                it.meta.ref_gff
                ]
        }
        .set { liftoff_in }

    LIFTOFF(liftoff_in, [])

    emit:
    ch_main
    quast_out               = QC.out.quast_out
    busco_out               = QC.out.busco_out
    merqury_report_files    = QC.out.merqury_report_files
}
