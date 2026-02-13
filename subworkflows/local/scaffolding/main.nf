include { RUN_LINKS         } from './links/main'
include { RUN_LONGSTITCH    } from './longstitch/main'
include { RUN_RAGTAG        } from './ragtag/main'
include { HIC               } from './hic/main'

workflow SCAFFOLD {
    take:
    ch_main
    meryl_kmers

    main:
    channel.empty().set { ch_versions }
    channel.empty().set { links_busco }
    channel.empty().set { links_quast }
    channel.empty().set { links_merqury }
    channel.empty().set { longstitch_busco }
    channel.empty().set { longstitch_quast }
    channel.empty().set { longstitch_merqury }
    channel.empty().set { ragtag_busco }
    channel.empty().set { ragtag_quast }
    channel.empty().set { ragtag_merqury }

    // There is no support for scaffolding of scaffolded scaffolds.
    // But it is possible that one sample is scaffolded with different tools.
    // Therefore main is filtered, instead of branched.

    ch_main
        .filter {
            it ->  it.meta.scaffold_links
        }
    .set { links_in }

    RUN_LINKS(links_in, meryl_kmers)
    RUN_LINKS.out.ch_main
        .set { links_out }

    ch_main
        .filter {
            it ->  it.meta.scaffold_longstitch
        }
    .set { longstitch_in }

    RUN_LONGSTITCH(longstitch_in, meryl_kmers)
    RUN_LONGSTITCH.out.ch_main
        .set { longstitch_out }

    ch_main
        .filter {
            it ->  it.meta.scaffold_hic
        }
    .set { hic_in }

    HIC(hic_in, meryl_kmers)
    HIC.out.ch_main
        .set { hic_out }

    ch_main
        .filter {
            it -> it.meta.scaffold_ragtag && !it.meta.hic_reads && !it.meta.scaffold_longstitch && !it.meta.scaffold_links
        }
    .mix(hic_out.filter { it -> it.meta.scaffold_ragtag } )
    .mix(longstitch_out.filter { it -> it.meta.scaffold_ragtag } )
    .mix(links_out.filter { it -> it.meta.scaffold_ragtag } )
    .set { ragtag_in }

    RUN_RAGTAG(ragtag_in, meryl_kmers)
    RUN_RAGTAG.out.ch_main
        .set { ragtag_out }

    // Deal with cases that are single scaffold
    links_out
        .filter {it -> !it.meta.scaffold_longstitch && !it.meta.scaffold_ragtag }
        .map { meta -> [ meta: meta - meta.subMap("links_scaffold") + [ scaffolds: [ links: meta.scaffolds_links ] ]  ]}
        .mix(
            longstitch_out
                .filter {it -> !it.meta.scaffold_links && !it.meta.scaffold_ragtag }
                .map { meta -> [ meta: meta - meta.subMap("scaffolds_longstitch") + [ scaffolds: [ longstitch: meta.scaffolds_longstitch ] ]  ]}
        )
        .mix(
            ragtag_out
                .filter {it -> !it.meta.scaffold_links && !it.meta.scaffold_longstitch }
                .map { meta -> [ meta: meta - meta.subMap("scaffolds_ragtag") + [ scaffolds: [ ragtag: meta.scaffolds_ragtag ] ]  ]}
        )
        .mix(
            hic_out
                .map { meta -> [ meta: meta - meta.subMap("scaffolds_hic") + [ scaffolds: [ hic: meta.scaffolds_hic ] ]  ]}

            )
        // mix in those that are double scaffolded: , links-ragtag, longstitch-ragtag
        // links-longstitch
        .mix(
            links_out
                .filter {it -> it.meta.scaffold_longstitch && !it.meta.scaffold_ragtag }
                .map {meta -> [meta.id, meta]}
                // Join without filtering, inner-join
                .join(
                    longstitch_out
                        .map {meta -> [meta.id, meta]}
                )
                .map {
                    _id, meta_links, meta_longstitch -> [
                        meta: meta_links -
                         meta_links.subMap("scaffolds_links") +
                         [scaffolds: [links: meta_links.scaffolds_links, longstitch: meta_longstitch.scaffolds_longstitch]] ]
                }
        )
        //links-ragtag
        .mix(
            links_out
                .filter {it -> !it.meta.scaffold_longstitch && it.meta.scaffold_ragtag }
                .map {meta -> [meta.id, meta]}
                // Join without filtering, inner-join
                .join(
                    ragtag_out
                        .map {meta -> [meta.id, meta]}
                )
                .map {
                    _id, meta_links, meta_ragtag -> [
                        meta: meta_links -
                         meta_links.subMap("scaffolds_links") +
                         [scaffolds: [links: meta_links.scaffolds_links, ragtag: meta_ragtag.scaffolds_ragtag]] ]
                }
        )
        //longstitch-ragtag
        .mix(
            longstitch_out
                .filter {it -> !it.meta.scaffold_links && it.meta.scaffold_ragtag }
                .map {meta -> [meta.id, meta]}
                // Join without filtering, inner-join
                .join(
                    ragtag_out
                        .map {meta -> [meta.id, meta]}
                )
                .map {
                    _id, meta_longstitch, meta_ragtag -> [
                        meta: meta_longstitch -
                         meta_longstitch.subMap("scaffolds_longstitch") +
                         [scaffolds: [longstitch: meta_longstitch.scaffolds_longstitch, ragtag: meta_ragtag.scaffolds_ragtag]] ]
                }
        )
        // mix in triple-scaffolded
        .mix(
            links_out
                .filter {it -> it.meta.scaffold_longstitch && it.meta.scaffold_ragtag }
                .map {meta -> [meta.id, meta]}
                // Join without filtering, inner-join
                .join(
                    longstitch_out
                        .map {meta -> [meta.id, meta]}
                )
                .join(
                    ragtag_out
                        .map {meta -> [meta.id, meta]}
                )
                .map {
                    _id, meta_links, meta_longstitch, meta_ragtag -> [
                        meta: meta_links -
                         meta_links.subMap("scaffolds_links") +
                            [
                                scaffolds: [
                                links: meta_links.scaffolds_links,
                                longstitch: meta_longstitch.scaffolds_longstitch,
                                ragtag: meta_ragtag.scaffolds_ragtag
                                ]
                            ]
                             ]
                }
        )
        .set { ch_main }


    RUN_LINKS.out.busco_out.set { links_busco }
    RUN_LINKS.out.quast_out.set { links_quast }
    RUN_LINKS.out.merqury_report_files.set { links_merqury }
    ch_versions = ch_versions.mix(RUN_LINKS.out.versions)

    RUN_LONGSTITCH.out.busco_out.set { longstitch_busco }
    RUN_LONGSTITCH.out.quast_out.set { longstitch_quast }
    RUN_LONGSTITCH.out.merqury_report_files.set { longstitch_merqury }
    ch_versions = ch_versions.mix(RUN_LONGSTITCH.out.versions)

    HIC.out.busco_out.set { hic_busco }
    HIC.out.quast_out.set { hic_quast }
    HIC.out.merqury_report_files.set { hic_merqury }
    ch_versions = ch_versions.mix(HIC.out.versions)

    RUN_RAGTAG.out.busco_out.set { ragtag_busco }
    RUN_RAGTAG.out.quast_out.set { ragtag_quast }
    RUN_RAGTAG.out.merqury_report_files.set { ragtag_merqury }
    ch_versions = ch_versions.mix(RUN_RAGTAG.out.versions)

    links_busco
        .concat(longstitch_busco)
        .concat(ragtag_busco)
        .concat(hic_busco)
        .set { scaffold_busco_reports }

    links_quast
        .concat(longstitch_quast)
        .concat(ragtag_quast)
        .concat(hic_quast)
        .set { scaffold_quast_reports }

    links_merqury
        .concat(longstitch_merqury)
        .concat(ragtag_merqury)
        .concat(hic_merqury)
        .set { scaffold_merqury_reports }

    versions = ch_versions

    emit:
    ch_main
    scaffold_busco_reports
    scaffold_quast_reports
    scaffold_merqury_reports
    versions
}
