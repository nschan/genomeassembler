include { COUNT } from '../../../../modules/local/jellyfish/count/main'
include { HISTO } from '../../../../modules/local/jellyfish/histo/main'
include { STATS } from '../../../../modules/local/jellyfish/stats/main'
include { GENOMESCOPE } from '../../../../modules/local/genomescope/main'

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
                    meta: [
                        id: it[1],
                        metas: it[0],
                        jellyfish_k: it[2][0],
                        qc_read_mean: it[4][0]
                    ],
                    qc_reads_path: it[3][0]
                ]
        }
        .mix(
            ch_main
                .filter { it -> !it.meta.group }
                .map {
                    it ->
                    [
                        meta: it.meta,
                        qc_reads_path: it.meta.qc_reads_path
                    ]
                }
        )

    COUNT(samples)
    kmers = COUNT.out.kmers

    HISTO(kmers)

    genomescope_in = HISTO.out.histo
        .map { meta, hist ->
                    [
                        meta,
                        hist,
                        meta.jellyfish_k,
                        meta.qc_read_mean

                    ]
        }

    STATS(kmers)

    GENOMESCOPE(genomescope_in)

    outputs = GENOMESCOPE.out.estimated_hap_len
        .filter { it -> it[0].metas }
        .flatMap { it ->
            it[0].metas
                .collect { meta -> [ meta: meta + [ genome_size: it[1] ] ] }
        }
        .mix(GENOMESCOPE.out.estimated_hap_len
            .filter { it -> !it[0].metas }
            .map {
                it -> [ meta: it[0] + [ genome_size: it[1] ] ]
            }
        )

    outputs.dump(tag: "Jellyfish outputs")

    genomescope_summary = GENOMESCOPE.out.summary

    genomescope_plot = GENOMESCOPE.out.plot

    emit:
    main_out = outputs
    genomescope_summary
    genomescope_plot
}
