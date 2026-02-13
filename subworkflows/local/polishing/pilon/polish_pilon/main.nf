include { RUN_PILON } from '../run_pilon/main'
include { MAP_SR } from '../../../mapping/map_sr/main'
include { RUN_LIFTOFF } from '../../../liftoff/main'
include { QC } from '../../../qc/main.nf'

workflow POLISH_PILON {
    take:
    ch_main
    meryl_kmers

    main:
    channel.empty().set { ch_versions }

    ch_main
        .multiMap {
            it ->
            shortreads: [it.meta, it.meta.shortreads]
            assembly: [
                it.meta,
                it.meta.polish == "medaka+pilon" ? it.meta.polished.medaka : it.meta.polish == "dorado+pilon" ? it.meta.polished.dorado : it.meta.assembly
                ]
        }
        .set { map_sr_in }

    //map_sr_in.shortreads.view {"POLISH_PILON: map_sr_in.shortreads: $it"}
    //map_sr_in.assembly.view {"POLISH_PILON: map_sr_in.assembly: $it"}

    MAP_SR(map_sr_in.shortreads, map_sr_in.assembly)

    RUN_PILON(map_sr_in.assembly, MAP_SR.out.aln_to_assembly_bam_bai)

    RUN_PILON.out.improved_assembly
        .set { pilon_polished }

    pilon_polished
        .map { meta, polished_pilon -> [ meta: meta + [ polished: [pilon: polished_pilon] ] ]  }
        .set { ch_main }

    ch_versions = ch_versions.mix(RUN_PILON.out.versions)

    QC(ch_main.map { it -> [meta: it.meta - it.meta.subMap("assembly_map_bam") + [assembly_map_bam: null] ]},
        pilon_polished.map {meta, polished -> [meta.id, polished ]},
        meryl_kmers)

    ch_versions = ch_versions.mix(QC.out.versions)

    ch_main
        .filter {
            it -> it.meta.lift_annotations
        }
        .map { it ->
                [
                it.meta,
                it.meta.polished.pilon,
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
