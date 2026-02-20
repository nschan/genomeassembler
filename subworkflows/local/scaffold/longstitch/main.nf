include { LONGSTITCH } from '../../../../modules/local/longstitch/main'
include { QC } from '../../qc/main'
include { LIFTOFF } from '../../../../modules/nf-core/liftoff/main'

workflow RUN_LONGSTITCH {
    take:
    ch_main
    meryl_kmers

    main:
    channel.empty().set { ch_versions }

    ch_main
        .map {
            it ->
            [
                it.meta,
                it.meta.polished ? (it.meta.polished.pilon ?: it.meta.polished.medaka ?: it.meta.polished.dorado) : it.meta.assembly,
                it.meta.qc_reads_path,
                it.meta.genome_size
            ]
        }
        .set { longstitch_in }

    longstitch_in.dump(tag: "SCAFFOLD: LONGSTITCH: inputs")

    LONGSTITCH(longstitch_in)

    LONGSTITCH.out.ntlLinks_arks_scaffolds
        .map { meta, scaff_longst -> [meta: meta + [scaffolds_longstitch: scaff_longst] ] }
        .set { ch_main_scaffolded }

    QC(ch_main_scaffolded.map { it -> [ meta: it.meta - it.meta.subMap("assembly_map_bam") + [assembly_map_bam: null] ] },
        LONGSTITCH.out.ntlLinks_arks_scaffolds.map { meta, scaffold -> [meta.id, scaffold]},
        meryl_kmers)

    ch_main_scaffolded
        .filter {
            it -> it.meta.lift_annotations
        }
        .map { it ->
                [
                it.meta,
                it.meta.scaffolds_longstitch,
                it.meta.ref_fasta,
                it.meta.ref_gff
                ]
        }
        .set { liftoff_in }

    LIFTOFF(liftoff_in,[])

    emit:
    ch_main                 = ch_main_scaffolded
    quast_out               = QC.out.quast_out
    busco_out               = QC.out.busco_out
    merqury_report_files    = QC.out.merqury_report_files
    versions                = ch_versions
}
