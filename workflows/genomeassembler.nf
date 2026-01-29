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
include { SCAFFOLD                  } from '../subworkflows/local/scaffolding/main'

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

    The "main" channel, contains all sample-wise information.
    This channel should be the main input of all subworkflows
    and the subworkflows should make changes to this map. The
    main channel should stay a map whenever possible and this
    main channel reflects all pipeline parameters.
    I will make use of the meta map to pass additional infor-
    mation into processes. This is neccessary to provide fine
    control for parameterization of processes. This is passed
    via ext.args to the process and fetched from meta.

    The keys are defined in
    ./subworkflows/local/utils_nfcore_genomeassembler/main.nf

        meta: [
            id: string,
            ontreads: path,
            hifireads: path,
            strategy: string,
            assembler_ont: string,
            assembler_hifi: string,
            scaffolding: string,
            genome_size: integer,
            assembler_ont_args: string,
            assembler_hifi_args: string,
            ref_fasta: path,
            ref_gff: path,
            shortread_F: path,
            shortread_R: path,
            paired: bool
            ont_collect: bool,
            ont_trim: bool,
            ont_jellyfish: bool,
            hifi_trim: bool,
            hifi_primers: path,
            polish_medaka: bool,
            medaka_model: string,
            polish_pilon: bool,
            scaffold_longstitch: bool,
            scaffold_links: bool,
            scaffold_ragtag: bool,
            use_ref: bool,
            flye_mode: string,
            assembly: path,
            ref_map_bam: path,
            assembly_map_bam: path,
            qc_reads: string ["ont","hifi"],
            qc_reads_path: path,
            quast: bool,
            busco: bool,
            busco_lineage: string,
            busco_db: path,
            lift_annotations: bool,
            shortread_F: path,
            shortread_R: path,
            paired: bool,
            use_short_reads: bool,
            shortread_trim: bool
        ]



    ===========
       JOINS
    ===========

    Since this channel needs to stay a map so I can pull out the correct elements, joining is difficult:
    Nextflow's join operator only works on list-typed channels, but the channels here are maps.
    For this reason, there are some confuding map operations involved where each map-element is converted to a list,
    containing the value and the previous map. The whole channel is turned into a list this way:

    Something like

    [
        meta: [id: something1],
        somepath: "/path"
    ]

    becomes

    [
        [id: something, meta: [id: something]],
        ["/path", somepath: "/path"]
    ]

    This can be joined to

    [
        [id: something, meta: [id: something]],
        ["different_path", otherpath: "different_path"]
    ]

    This makes it possible to join on the first element (the one containing meta):

    [
        [id: something, meta: [id: something]],
        ["path", somepath: "path"],
        ["different_path", otherpath: "different_path"]
    ]

    After joining, the map is recreated from the second list element, to create:

    [
        meta: [id: something],
        somepath: "path",
        otherpath: "different_path"
    ]

    This is (sadly) a somewhat frequent pattern in this pipeline and
    it is done like this:

    map_channel_1
            // Convert to list for join
            .map { it -> it.collect { entry -> [ entry.value, entry ] } }
            .join( map_channel_2
                     // Convert to list for join
                    .map { it -> it.collect { entry -> [ entry.value, entry ] } }
            )
            // After joining re-create the maps from the stored map
            .map { it -> it.collect { _entry, map -> [ (map.key): map.value ] }.collectEntries() }
    */
    channel.empty().set { meryl_kmers }

    channel.empty().set { ch_versions }

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


    //ch_main_prepared.view{"Main WF: Prepared out: $it "}
    /*
    Assembly
    */
    // This pipeline is named genomeassembler, so everything goes into assemble
    // even it might not actually be assembled.
    ASSEMBLE(ch_main_prepared, meryl_kmers)

    ASSEMBLE.out.ch_main.set { ch_main_assembled }

    ch_versions = ch_versions.mix(ASSEMBLE.out.versions)
    /*
    Polishing
    */
    ch_main_assembled
        .branch {
            it ->
            def polishers = ["pilon", "medaka", "medaka+pilon", "dorado", "dorado+pilon"]
            /*debug
            def polishValue = it.meta.polish
            def inList = polishers.contains(polishValue)
            println "DEBUG: polish='${polishValue}' (type: ${polishValue.class.name}), inList=${inList}"
            polish:     inList
            no_polish:  !inList
            DEBUG: polish='"medaka+pilon"' (type: java.lang.String), inList=false
            No quotes in samplesheet?
            */
            polish:     polishers.contains(it.meta.polish) == true
            no_polish:  true
        }
        .set { ch_main_assembled_branched }

    ch_main_assembled_branched.polish.view {"ch_main_assembled_branched.polish: $it"}
    ch_main_assembled_branched.no_polish.view {"ch_main_assembled_branched.no_polish: $it"}

    POLISH(ch_main_assembled_branched.polish, meryl_kmers)

    ch_main_assembled_branched.no_polish
        .mix(POLISH.out.ch_main)
        .set { ch_main_polished }

    ch_versions = ch_versions.mix(POLISH.out.versions)

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

    ch_versions = ch_versions.mix(PREPARE.out.versions).mix(ASSEMBLE.out.versions).mix(POLISH.out.versions).mix(SCAFFOLD.out.versions)

    ch_versions = ch_versions



    /*
    Report
    */
    softwareVersionsToYAML(ch_versions)
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

    //fastplong_jsons.view { it -> "UNQIE JSONS: $it"}

    REPORT( report_files,
            report_functions,
            report_scripts,
            fastplong_jsons,
            genomescope_files,
            quast_files,
            busco_files,
            merqury_files,
            channel.fromPath("${params.outdir}/pipeline_info/nf_core_pipeline_software_versions.yml"),
            ch_main.map { it -> [sample: [id: it.meta.id, group: it.group]]}.collect()
    )

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
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

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'genomeassembler_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    _report = REPORT.out.report_html.toList()

    emit:
    _report
    versions = ch_versions // channel: [ path(versions.yml) ]
}
