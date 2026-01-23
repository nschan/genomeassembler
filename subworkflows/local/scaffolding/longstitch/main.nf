include { LONGSTITCH } from '../../../../modules/local/longstitch/main'
include { QC } from '../../qc/main'
include { RUN_LIFTOFF } from '../../liftoff/main'

workflow RUN_LONGSTITCH {
    take:
    ch_main
    meryl_kmers

    /*
    TODO:
    Longstitch needs genomesize. For ONT reads that is estimated if not provided.
    For hifireads it needs to be provided. Currently, not checks for that..
    Depending on the reads used, longmap needs to be changed.
    This should probably be done via args / config, but needs a way to do this per-sample.
    Probably passing in additional information via meta is the way to go for this.
    */

    main:
    channel.empty().set { ch_versions }

    ch_main
        .map {
            it ->
            [
                it.meta,
                it.meta.polished ? (it.polished.pilon ?: it.polished.medaka) : it.assembly,
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

    ch_versions = ch_versions.mix(LONGSTITCH.out.versions)

    QC(ch_main_scaffolded.map { it -> [ meta: it.meta - it.meta.subMap("assembly_map_bam") + [assembly_map_bam: null] ] },
        LONGSTITCH.out.ntlLinks_arks_scaffolds.map { meta, scaffold -> [meta.id, scaffold]},
        meryl_kmers)

    ch_versions = ch_versions.mix(QC.out.versions)

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

    RUN_LIFTOFF(liftoff_in)
    ch_versions = ch_versions.mix(RUN_LIFTOFF.out.versions)

    emit:
    ch_main
    quast_out               = QC.out.quast_out
    busco_out               = QC.out.busco_out
    merqury_report_files    = QC.out.merqury_report_files
    versions                = ch_versions
}
