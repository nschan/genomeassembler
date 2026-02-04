include { FASTP                 } from '../../../../modules/nf-core/fastp/main'
include { FASTP as FASTP_HIC    } from '../../../../modules/nf-core/fastp/main'
include { MERYL_COUNT           } from '../../../../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM        } from '../../../../modules/nf-core/meryl/unionsum/main'

workflow PREPARE_SHORTREADS {
    take:
    shortreads_in

    main:
    channel.empty().set { ch_versions }

    shortreads_in
        .map { row -> row.meta.shortread_F ? create_shortread_channel(row.meta) : row } // function below
        .branch {
            it ->
                trim: it.meta.shortread_trim
                no_trim: !it.meta.shortread_trim
        }
        .set { shortreads }

    shortreads_in
        .map { row -> row.meta.hic_F ? create_hic_shortread_channel(row.meta) : row }
        .branch {
            row ->
                trim: row.meta.hic_trim
                no_trim: !row.meta.hic_trim
        }
        .set { hic_trim }

    shortreads.trim.dump(tag: "shortread trim channel")
    hic_trim.trim.dump(tag: "hic trim channel")

    shortreads
        .trim
        .filter { it -> it.meta.group }
        .map { it -> [it.meta, it.meta.group] }
        .groupTuple(by: 1)
        .map {
            it ->
                [
                    [
                        id: it[1], // the group
                        metas: it[0]
                    ],
                    it[0].shortreads[0], // Pull path from meta
                    []
                ]
        }
        .mix(
            shortreads
                .trim
                .filter { it -> !it.meta.group }
                .map {
                it -> [ it.meta, it.meta.shortreads, [] ]
                }
        )
        .set { trim_in }

    hic_trim
        .trim
        .filter { it -> it.meta.group }
        .map {it -> [it.meta, it.meta.group]}
        .groupTuple(by: 1)
        .map {
            it ->
                [
                    [
                        id: it[1], // the group
                        metas: it[0]
                    ],
                    it[0].hic_reads[0], // Pull path from meta
                    []
                ]
        }
        .mix(
            hic_trim
                .trim
                .filter { it -> !it.meta.group }
                .map {
                    it -> [ it.meta, it.meta.hic_reads, [] ]
                }
        )
        .set { hic_trim_in }

    trim_in.dump(tag: "Trim in")

    FASTP(trim_in, false, false, false)

    FASTP.out.reads
        .filter { it -> it[0].metas }
        .flatMap { it -> // looks like [meta <[id, metas]>, output_path]
            it[0].metas
                  .collect { meta -> [ meta: meta - meta.subMap("shortreads") + [ shortreads: it[1] ] ] }
        }
        .mix(
            FASTP.out.reads
                .filter { it -> !it[0].metas }
                .map { it -> [ meta: it[0] - it[0].subMap("shortreads") + [ shortreads: it[1] ] ] }
        )
        .set { trimmed_reads }

    trimmed_reads.dump(tag: "Trim out")
    // unite branched:
    // add trimmed reads to trim channel, then mix with shortreads.no_trim

    FASTP_HIC(hic_trim_in, false, false, false)

    FASTP_HIC.out.reads
        .filter { it -> it[0].metas }
        .flatMap { it -> // looks like [meta <[id, metas]>, output_path]
            it[0].metas
                  .collect { meta -> [ meta: meta - meta.subMap("hic_reads") + [ hic_reads: it[1] ] ] }
        }
        .mix(
            FASTP_HIC.out.reads
                .filter { it -> !it[0].metas }
                .map { it -> [ meta: it[0] - it[0].subMap("hic_reads") + [ hic_reads: it[1] ] ] }
        )
        .set { hic_trimmed_reads }

    trimmed_reads
        .mix( shortreads.no_trim )
        .set { shortreads }

    // add HiC trimmed to those that need it

    shortreads
        .filter { row -> row.meta.hic_trim }
        .map { row -> [ row.meta.id, row.meta ] }
        .join(
            hic_trimmed_reads
                .map { meta  ->
                    [
                        meta.id,
                        meta.hic_reads
                    ]
                }
        )
        .map {
            _id, meta, trimmed_hic_reads ->
                [
                    meta: meta - meta.subMap("hic_reads") + [ hic_reads: trimmed_hic_reads ]
                ]
        }
        .mix(
            trimmed_reads
                .filter { row -> !row.meta.hic_trim }
        )
        .set { shortreads }


    ch_versions = ch_versions.mix(FASTP.out.versions)

    shortreads
        .filter { it -> it.meta.merqury }
        .filter { it -> it.meta.group  }
        .map { it -> [ it.meta, it.meta.group, it.meta.shortreads, it.meta.meryl_k ] }
        // Create a group
        .groupTuple(by: 1)
        .map {
            it -> [
                meta: [ id: it[1], metas: it[0] ],
                shortreads: it[2][0],
                meryl_k: it[3][0]
            ]
        }
        .mix(shortreads
            .filter { it -> it.meta.merqury }
            .filter { it -> !it.meta.group  }
            .map { it -> [meta: it.meta, shortreads: it.meta.shortreads, meryl_k: it.meta.meryl_k]}
        )
        .multiMap { it ->
            reads: [ it.meta, it.shortreads ]
            kmer_size: it.meryl_k
        }
        .set { meryl_in }

    MERYL_COUNT(meryl_in.reads, meryl_in.kmer_size)

    MERYL_UNIONSUM(MERYL_COUNT.out.meryl_db, params.meryl_k)

    MERYL_UNIONSUM.out.meryl_db
        .filter { it -> it[0].metas }
        .flatMap { it -> // looks like [meta <[id, metas]>, output_path]
            it[0].metas
                  .collect { meta -> [ meta, it[1] ] }
        }
        .mix(MERYL_UNIONSUM.out.meryl_db
            .filter { it -> !it[0].metas }
            .map {
                it -> [ it[0], it[1] ]
            }
        )
        .map {meta , kmers -> [meta.id, kmers]}
        .set { meryl_kmers }



    versions = ch_versions.mix(MERYL_COUNT.out.versions).mix(MERYL_UNIONSUM.out.versions)

    emit:
    main_out        = shortreads
    meryl_kmers
    versions
    fastp_json      = FASTP.out.json
}

def create_shortread_channel(row) { // This function expects a meta map as input
    // create meta map
    def meta = row
    meta.paired = row.paired.toBoolean()
    meta.single_end = !meta.paired

    // add path(s) of the fastq file(s) to the meta map
    def shortreads = []
    if (!file(row.shortread_F).exists()) {
        exit(1, "ERROR: shortread_F fastq file does not exist!\n${row.shortread_F}")
    }
    if (!meta.paired) {
        shortreads = [meta: meta + [shortreads: [file(row.shortread_F)]]]
    }
    else {
        if (!file(row.shortread_R).exists()) {
            exit(1, "ERROR: shortread_R fastq file does not exist!\n${row.shortread_R}")
        }
        shortreads = [ meta: meta + [shortreads:  [file(row.shortread_F), file(row.shortread_R)]] ]
    }
    return shortreads
}

def create_hic_shortread_channel(row) { // This function expects a meta map as input
    // create meta map
    def meta = row
    meta.paired = true
    meta.single_end = !meta.paired

    // add path(s) of the fastq file(s) to the meta map
    def hic_reads = []
    if (!file(row.hic_F).exists()) {
        exit(1, "ERROR: hic_F fastq file does not exist!\n${row.hic_F}")
    }
    if (!meta.paired) {
        hic_reads = [meta: meta + [hic_reads: [file(row.hic_F)]]]
    }
    else {
        if (!file(row.hic_R).exists()) {
            exit(1, "ERROR: shortread_R fastq file does not exist!\n${row.hic_R}")
        }
        hic_reads = [ meta: meta + [hic_reads:  [file(row.hic_F), file(row.hic_R)]] ]
    }
    return hic_reads
}
