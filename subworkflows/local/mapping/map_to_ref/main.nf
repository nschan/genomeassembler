include { MINIMAP2_ALIGN as ALIGN } from '../../../../modules/nf-core/minimap2/align/main'
include { BAM_STATS_SAMTOOLS as BAM_STATS } from '../../../nf-core/bam_stats_samtools/main'

workflow MAP_TO_REF {
    take:
    ch_map_ref // meta, reads, refs

    main:
    // Map reads to reference
    ALIGN(ch_map_ref, true, 'bai', false, false)

    ALIGN.out.bam
        .set { ch_aln_to_ref_bam }

    ALIGN.out.index
        .set { aln_to_ref_bai }

    ch_aln_to_ref_bam
        .join(aln_to_ref_bai)
        .set { ch_aln_to_ref_bam_bai }

    ch_map_ref
        .map { meta, _reads, fasta -> [[meta], fasta] }
        .set { ch_fasta }

    BAM_STATS(ch_aln_to_ref_bam_bai, ch_fasta)

    emit:
    ch_aln_to_ref_bam //  meta, bam
}
