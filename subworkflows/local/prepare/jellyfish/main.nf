include { JELLYFISH_COUNT as COUNT } from '../../../../modules/nf-core/jellyfish/count/main'
include { HISTO } from '../../../../modules/local/jellyfish/histo/main'
include { STATS } from '../../../../modules/local/jellyfish/stats/main'
include { GENOMESCOPE2 } from '../../../../modules/nf-core/genomescope2/main'

workflow JELLYFISH {
    take:
    ch_main

    main:
    genomescope_in = channel.empty()

    samples = ch_main
        .filter { it -> it.meta.group }
        .map { it ->
            [
                it.meta,
                it.meta.group,
                it.meta.jellyfish_k,
                it.meta.qc_reads_path,
                it.meta.qc_read_mean
            ]
        }
        .groupTuple(by: 1)
        .map {
            it ->
                [
                    [
                        id: it[1],
                        metas: it[0],
                        jellyfish_k: it[2][0],
                        qc_read_mean: it[4][0]
                    ],
                    it[3][0]
                ]
        }
        .mix(
            ch_main
                .filter { it -> !it.meta.group }
                .map {
                    it ->
                    [
                    it.meta,
                    it.meta.qc_reads_path
                    ]
                }
        )
    jellyfish_count_in = samples
        .multiMap { meta, reads ->
            fasta: [meta, reads]
            kmer_length: meta.jellyfish_k
            size:  meta.jellyfish_size
            }

    COUNT(jellyfish_count_in.fasta, jellyfish_count_in.kmer_length, jellyfish_count_in.size)
    kmers = COUNT.out.jf

    HISTO(kmers)

    genomescope_in = HISTO.out.histo

    STATS(kmers)

    GENOMESCOPE2(genomescope_in)

    outputs = GENOMESCOPE2.out.estimated_hap_len
        .filter { it -> it[0].metas }
        .flatMap { it ->
            it[0].metas
                .collect { meta -> [ meta: meta + [ genome_size: it[1] ] ] }
        }
        .mix(GENOMESCOPE2.out.estimated_hap_len
            .filter { it -> !it[0].metas }
            .map {
                it -> [ meta: it[0] + [ genome_size: it[1] ] ]
            }
        )

    outputs.dump(tag: "Jellyfish outputs")

    genomescope_summary = GENOMESCOPE2.out.summary

    genomescope_plot = GENOMESCOPE2.out.transformed_log_plot_png

    emit:
    main_out = outputs
    genomescope_summary
    genomescope_plot
}
