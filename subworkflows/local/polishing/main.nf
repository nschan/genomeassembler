include { POLISH_MEDAKA } from './medaka/polish_medaka/main.nf'
include { POLISH_PILON } from './pilon/polish_pilon/main.nf'
include { POLISH_DORADO } from './dorado/main.nf'

workflow POLISH {
    take:
    ch_main
    meryl_kmers

    main:
    channel.empty().set { polish_busco_reports }
    channel.empty().set { polish_quast_reports }
    channel.empty().set { polish_merqury_reports }

    ch_main
        .branch { it ->
            def medaka_polishers = ["medaka","medaka+pilon"]
            def dorado_polishers = ["dorado","dorado+pilon"]
            medaka: medaka_polishers.contains(it.meta.polish)
            dorado: dorado_polishers.contains(it.meta.polish)
            no_ont_polish: !medaka_polishers.contains(it.meta.polish) && !dorado_polishers.contains(it.meta.polish)
        }
        .set { ch_main_polish }

    POLISH_MEDAKA(ch_main_polish.medaka, meryl_kmers)

    POLISH_DORADO(ch_main_polish.dorado, meryl_kmers)

    POLISH_MEDAKA.out.busco_out
        .mix(POLISH_DORADO.out.busco_out)
        .set { polish_busco_reports }

    POLISH_MEDAKA.out.quast_out
        .mix(POLISH_DORADO.out.quast_out)
        .set { polish_quast_reports }

    POLISH_MEDAKA.out.merqury_report_files
        .mix(POLISH_DORADO.out.merqury_report_files)
        .set { polish_merqury_reports }

    POLISH_MEDAKA.out.ch_main
        .mix(POLISH_DORADO.out.ch_main)
        .mix(ch_main_polish.no_ont_polish)
        .set { ch_main_polish_pilon }



    /*
    Polishing with short reads using pilon
    */

    ch_main_polish_pilon
        .branch {
            it ->
            def pilon_polishers = ["pilon","medaka+pilon", "dorado+pilon"]
            pilon: pilon_polishers.contains(it.meta.polish)
            no_pilon: true
        }
        .set { ch_main_polish_pilon_in }

    POLISH_PILON(ch_main_polish_pilon_in.pilon, meryl_kmers)

    ch_main_polish_pilon_in.no_pilon.mix(POLISH_PILON.out.ch_main)
        .set { ch_out }

    polish_busco_reports
        .concat(
            POLISH_PILON.out.busco_out
        )
        .set { polish_busco_reports }

    polish_quast_reports
        .concat(
            POLISH_PILON.out.quast_out
        )
        .set { polish_quast_reports }

    polish_merqury_reports
        .concat(
            POLISH_PILON.out.merqury_report_files
        )
        .set { polish_merqury_reports }

    emit:
    ch_main = ch_out
    polish_busco_reports
    polish_quast_reports
    polish_merqury_reports
}
