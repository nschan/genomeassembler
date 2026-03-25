include { MAP_TO_ASSEMBLY } from '../mapping/map_to_assembly/main'
include { QUAST } from '../../../modules/local/quast/main'
include { BUSCO_BUSCO as BUSCO } from '../../../modules/nf-core/busco/busco/main'
include { MERQURY_MERQURY as MERQURY } from '../../../modules/nf-core/merqury/merqury/main'

workflow QC {
    take:
    ch_main
    scaffolds
    meryl_kmers

    main:
    channel.empty().set { quast_out }
    channel.empty().set { busco_out }
    channel.empty().set { merqury_report_files }

    ch_main
        .branch {
            it ->
            shortread: it.meta.use_short_reads
            no_shortread: !it.meta.use_short_reads
        }
        .set { ch_shortread_branched }

    ch_shortread_branched
        .shortread
        .filter { it -> it.meta.merqury }
        .map { it -> [it.meta.id, it.meta] }
        .join(scaffolds)
        .join(meryl_kmers)
        .map { _id, meta, scaffs, kmers ->
                [ meta, kmers, scaffs ]
            }
        .set { merqury_in }

    MERQURY(merqury_in)

    // Make sure that Polish and Scaffold main channels do not contain assembly_map_bam

    ch_main
        .branch {
            it ->
            map_to_assembly: it.meta.quast && !it.meta.assembly_map_bam
            no_map_to_assembly: !it.meta.quast || (it.meta.quast && it.meta.assembly_map_bam)
        }
        .set { ch_map_branched }

    ch_map_branched
        .map_to_assembly
        .map {
            it -> [ it.meta.id, it.meta ]
        }
        .join(scaffolds)
        .map {
            _id, meta, target_scaffolds ->
            [
                meta + [qc_target: target_scaffolds], // QC Target only exists in QC channel, and takes the scaffold that should be qc'ed
                meta.qc_reads_path,
                target_scaffolds
            ]
        }
        .set { map_assembly_in }

    MAP_TO_ASSEMBLY(map_assembly_in)

     // create main channel with mappings
    MAP_TO_ASSEMBLY.out.aln_to_assembly_bam
        .map { meta, assembly_map_bam ->
            [
                meta: meta + [ assembly_map_bam: assembly_map_bam ]
            ]
        }
        .mix(ch_map_branched.no_map_to_assembly)
        .set { ch_qc }

    ch_qc
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
                use_gff: it.meta.use_ref && it.meta.ref_gff ? true : false
            }
        .set { quast_in }

    QUAST(quast_in.quast_in, quast_in.use_ref, quast_in.use_gfgf)
    QUAST.out.tsv.set { quast_out }

    ch_qc
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
    BUSCO.out.batch_summary.set { busco_out }


    MERQURY.out.stats
        .join(
            MERQURY.out.spectra_asm_hist
        )
        .join(
            MERQURY.out.spectra_cn_hist
        )
        .join(
            MERQURY.out.assembly_qv
        )
        .set { merqury_report_files }

    emit:
    ch_main     // QC does not (and should not) modify ch_main but returns the input.
    quast_out
    busco_out
    merqury_report_files
}
