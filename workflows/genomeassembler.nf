/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_genomeassembler_pipeline'

// Read preparation
include { PREPARE                   } from '../subworkflows/local/prepare/main'

// Assembly
include { ASSEMBLE                  } from '../subworkflows/local/assemble/main'

// Polishing
include { POLISH                    } from '../subworkflows/local/polishing/main'

// Scaffolding
include { SCAFFOLD                  } from '../subworkflows/local/scaffold/main'

// reporting
include { REPORT                    } from '../modules/local/report/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMEASSEMBLER {
    take:
    ch_input

    main:
    // Initialize empty channels
    ch_input.set { ch_main }

    /*
    This pipeline uses a "meta-stuffing" appraoch. All information
    about a sample is always stored in a map stored in [0]/"meta".
    Values are extracted from the map to create input channels.
    The correspoding key is created or updated from outputs.
    This largely eliminates the need for joins.

    The initial keys are defined in
    ./subworkflows/local/utils_nfcore_genomeassembler/main.nf
    */
    channel.empty().set { meryl_kmers }

    // Initialize channels for QC report collection
    channel
        .of([])
        .tap { quast_files }
        .tap { fastplong_jsons }
        .tap { genomescope_files }
        .map { it -> ["dummy", it] }
        .tap { busco_files }
        .map { it -> [it[0], it[1], it[1], it[1], it[1]] }
        .tap { merqury_files }

    /*
    =============
    Prepare reads
    =============
    */
    PREPARE(ch_main)

    PREPARE.out.ch_main.set { ch_main_prepared }

    PREPARE.out.meryl_kmers.set { meryl_kmers }

    /*
    Assembly
    */
    // This pipeline is named genomeassembler, so everything goes into assemble
    // even it might not actually be assembled.

    ASSEMBLE(ch_main_prepared, meryl_kmers)

    ASSEMBLE.out.ch_main.set { ch_main_assembled }

    /*
    Polishing
    */
    ch_main_assembled
        .branch {
            it ->
            def polishers = ["pilon", "medaka", "medaka+pilon", "dorado", "dorado+pilon"]
            polish:     polishers.contains(it.meta.polish)
            no_polish:  true
        }
        .set { ch_main_assembled_branched }

    POLISH(ch_main_assembled_branched.polish, meryl_kmers)

    ch_main_assembled_branched.no_polish
        .mix(POLISH.out.ch_main)
        .set { ch_main_polished }
    // Update scaffold for meta map

    ch_main_polished
        .branch { it ->
            scaffold: it.meta.scaffold_links || it.meta.scaffold_longstitch || it.meta.scaffold_ragtag
            no_scaffold: !it.meta.scaffold_links && !it.meta.scaffold_longstitch && !it.meta.scaffold_ragtag
        }
    .set {
        ch_main_polished_branched
    }
    /*
    Scaffolding
    */
    SCAFFOLD(ch_main_polished_branched.scaffold, meryl_kmers)

    // Recreate ch_main, even though it is not used since there are no later steps.

    ch_main_polished_branched
        .no_scaffold
        .mix(SCAFFOLD.out.ch_main)
        .set { ch_main_scaffolded }

    PREPARE.out.fastplong_json_reports
        .map { it -> it[1] }
        .unique()
        .collect()
        .set { fastplong_jsons }

    PREPARE.out.genomescope_summary
        .concat(
            PREPARE.out.genomescope_plot
        )
        .unique()
        .collect { it -> it[1] }
        .set { genomescope_files }

    def topic_versions = channel.topic("versions")
      .distinct()
      .branch { entry ->
          versions_file: entry instanceof Path
          versions_tuple: true
      }

    def topic_versions_string = topic_versions.versions_tuple
      .map { process, tool, version ->
          [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
      }
      .groupTuple(by:0)
      .map { process, tool_versions ->
          tool_versions.unique().sort()
          "${process}:\n${tool_versions.join('\n')}"
      }
    ch_collated_versions = topic_versions_string
    /*
    Report
    */
    ch_collated_versions
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'pipeline_software_' + 'versions.yml',
            sort: true,
            newLine: true
        )

    quast_files
        .mix(
            ASSEMBLE.out.assembly_quast_reports
            .mix(
                POLISH.out.polish_quast_reports
            )
            .mix(
                SCAFFOLD.out.scaffold_quast_reports
            )
        )
        .unique()
        .collect()
        .set { quast_files }

    busco_files
        .mix(
            ASSEMBLE.out.assembly_busco_reports
            .mix(
                POLISH.out.polish_busco_reports
            )
            .mix(
                SCAFFOLD.out.scaffold_busco_reports
            )
        )
        .unique()
        .collect { it -> it[1] }
        .set { busco_files }

    merqury_files
        .mix(
            ASSEMBLE.out.assembly_merqury_reports
            .mix(
                POLISH.out.polish_merqury_reports
            )
            .mix(
                SCAFFOLD.out.scaffold_merqury_reports
            )
        )
        .collect { it -> [it[1], it[2], it[3], it[4]] }
        .toSet()
        .flatten()
        .collect()
        .set { merqury_files }

    channel
        .fromPath("${projectDir}/assets/report/*")
        .collect()
        .set { report_files }
    // Report files
    channel
        .fromPath("${projectDir}/assets/report/functions/*")
        .collect()
        .set { report_functions }
    channel
        .fromPath("${projectDir}/assets/report/scripts/*")
        .collect()
        .set { report_scripts }

    REPORT( report_files,
            report_functions,
            report_scripts,
            fastplong_jsons,
            genomescope_files,
            quast_files,
            busco_files,
            merqury_files,
            channel.fromPath("${params.outdir}/pipeline_info/nf_core_pipeline_software_versions.yml"),
            ch_main.map { it -> [sample: [id: it.meta.id, group: it.meta.group]] }.collect()
    )

    _report = REPORT.out.report_html.toList()

    emit:
    _report
}
