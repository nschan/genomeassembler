include { COUNT } from '../../../../modules/local/jellyfish/count/main'
include { HISTO } from '../../../../modules/local/jellyfish/histo/main'
include { STATS } from '../../../../modules/local/jellyfish/stats/main'
include { GENOMESCOPE } from '../../../../modules/local/genomescope/main'

workflow JELLYFISH {
    take:
    ch_main

    main:
    channel.empty().set { genomescope_in }
    channel.empty().set { ch_versions }

    ch_main
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
                        jellyfish_k: it[2].unique()[0],
                        qc_read_mean: it[4].unique()[0]
                    ],
                    qc_reads_path: it[3].unique()[0]
                ]
        }
        .mix(
            ch_main
                .filter { it -> !it.group }
                .map {
                    it ->
                    [
                        meta: it.meta,
                        qc_reads_path: it.meta.qc_reads_path
                    ]
                }
        )
        .set { samples }

    COUNT(samples)
    COUNT.out.kmers.set { kmers }

    ch_versions = ch_versions.mix(COUNT.out.versions)

    HISTO(kmers)
    ch_versions = ch_versions.mix(HISTO.out.versions)

    HISTO.out.histo
        .map { meta, hist ->
                    [
                        meta,
                        meta.jellyfish_k,
                        meta.qc_read_mean,
                        hist
                    ]
        }
        .set { genomescope_in }

    STATS(kmers)

    GENOMESCOPE(genomescope_in)

    ch_versions = ch_versions
        .mix(GENOMESCOPE.out.versions)
        .mix(STATS.out.versions)

    GENOMESCOPE.out.estimated_hap_len
        .filter { it -> it[0].metas }
        .flatMap { it ->
            it[0].metas
                .collect { meta -> [ meta: meta + [ genome_size: it[1] ] ] }
        }
        .mix(GENOMESCOPE.out.estimated_hap_len
            .filter { it -> !it[0].ids }
            .map {
                it -> [ meta: it[0] + [ genome_size: it[1] ] ]
            }
        )
        .set { outputs }

    outputs.dump(tag: "Jellyfish outputs")

    GENOMESCOPE.out.summary.set { genomescope_summary }

    GENOMESCOPE.out.plot.set { genomescope_plot }

    versions = ch_versions

    emit:
    main_out = outputs
    genomescope_summary
    genomescope_plot
    versions
}
