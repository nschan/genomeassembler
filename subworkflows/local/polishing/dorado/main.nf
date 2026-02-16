include { DORADO_ALIGNER as ALIGN } from '../../../../modules/local/dorado/aligner/main.nf'
include { DORADO_POLISH as POLISH } from '../../../../modules/local/dorado/polish/main.nf'
include { QC } from '../../qc/main.nf'
include { LIFTOFF } from '../../../../modules/nf-core/liftoff/main'

workflow POLISH_DORADO {
    take:
    ch_main
    meryl_kmers

    main:

    ch_main
        .map { it -> [it.meta, it.meta.assembly, it.meta.ontreads] }
        .set { ch_aln_in }

    ALIGN(ch_aln_in)

    ALIGN.out.bam
        .join(ALIGN.out.bai)
        .map {meta, bam, bai-> [ meta, meta.assembly, bam, bai ] }
        .set { ch_polish_in }

    POLISH(ch_polish_in, [])

    POLISH.out.polished_alignment.set { polished_assembly }

    polished_assembly
        .map { meta, polished_dorado -> [meta: meta + [ polished: [polished_dorado: polished_dorado ] ] ]}
        .set { ch_main_out }

    QC(
        ch_main_out.map { it -> [meta: it.meta - it.meta.subMap("assembly_map_bam") + [ assembly_map_bam: null] ] },
        polished_assembly.map { meta, polished -> [meta.id, polished] },
        meryl_kmers
    )

    ch_main_out
        .filter {
            it -> it.meta.lift_annotations
        }
        .map { it ->
                [
                it.meta,
                it.meta.polished.dorado,
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
