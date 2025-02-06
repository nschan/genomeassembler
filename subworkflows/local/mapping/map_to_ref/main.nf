include { MINIMAP2_ALIGN as ALIGN } from '../../../../modules/nf-core/minimap2/align/main'
include { BAM_INDEX_STATS_SAMTOOLS as BAM_STATS } from '../../bam_sort_stat/main'

workflow MAP_TO_REF {
    take:
    in_reads
    ch_refs

    main:
    // Map reads to reference
    in_reads
        .join(ch_refs)
        .set { ch_map_ref_in }

    ALIGN(ch_map_ref_in ,true, false, false, false)

    ALIGN.out.bam.set { ch_aln_to_ref }
    ch_map_ref_in
        .map { meta, _reads, fasta -> [ meta, fasta ] }
        .set { ch_fasta }

    BAM_STATS(ch_aln_to_ref, ch_fasta)

    emit:
    ch_aln_to_ref
}
