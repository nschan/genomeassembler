include { LINKS } from '../../../../modules/nf-core/links/main'
include { QC } from '../../qc/main'
include { RUN_LIFTOFF } from '../../liftoff/main'

workflow RUN_LINKS {
    take:
    ch_main
    meryl_kmers

    main:
    channel.empty().set { ch_versions }
    ch_main.dump(tag: "SCAFFOLD: LINKS: WORKFLOW inputs")
    ch_main
        .multiMap { it ->
            assembly:   [it.meta, it.meta.polished ? (it.meta.polished.pilon ?: it.meta.polished.medaka ?: it.meta.polished.dorado) : it.meta.assembly]
            reads:      [it.meta, it.meta.qc_reads_path]
        }
        .set { links_in }

    links_in.assembly.dump(tag: "SCAFFOLD: LINKS: Assembly inputs")
    links_in.reads.dump(tag: "SCAFFOLD: LINKS: Read inputs")

    LINKS(links_in.assembly, links_in.reads)
    LINKS.out.scaffolds_fasta
        .map { meta, scaff_links -> [meta: meta + [scaffolds_links: scaff_links] ] }
        .set { ch_main_scaffolded }

    ch_versions = ch_versions.mix(LINKS.out.versions)

    QC(ch_main_scaffolded.map { it -> [meta: it.meta - it.meta.subMap("assembly_map_bam") + [assembly_map_bam: null] ]},
        LINKS.out.scaffolds_fasta.map { meta, scaffold -> [meta.id, scaffold]},
         meryl_kmers)

    ch_versions = ch_versions.mix(QC.out.versions)

    ch_main_scaffolded
        .filter {
            it -> it.lift_annotations
        }
        .map { it ->
                [
                it.meta,
                it.scaffolds_links,
                it.ref_fasta,
                it.ref_gff
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
