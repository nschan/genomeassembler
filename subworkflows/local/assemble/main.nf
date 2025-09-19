include { FLYE } from '../../../modules/nf-core/flye/main'
include { HIFIASM } from '../../../modules/nf-core/hifiasm/main'
include { HIFIASM as HIFIASM_ONT } from '../../../modules/nf-core/hifiasm/main'
include { GFA_2_FA as GFA_2_FA_HIFI } from '../../../modules/local/gfa2fa/main'
include { GFA_2_FA as GFA_2_FA_ONT} from '../../../modules/local/gfa2fa/main'
include { MAP_TO_REF } from '../mapping/map_to_ref/main'
include { RUN_LIFTOFF } from '../liftoff/main'
include { RAGTAG_PATCH } from '../../../modules/nf-core/ragtag/patch/main'
include { QC } from '../qc/main'


workflow ASSEMBLE {
    take:
    ch_main
    meryl_kmers

    main:
    // Empty channels
    Channel.empty().set { ch_versions }

    /*
    Samples are split into those that need assembly, and those that will not be assembled (i.e. assemblies are provided)
    */
    ch_main.dump(tag: "Assemble - Inputs")
    ch_main
        .branch {
            it ->
            to_assemble: !it.assembly
            no_assemble: it.assembly
        }
        .set {
            ch_main_branched
        }
    /*
    There are three assembly strategies:
        - Single: Using a single assembler with one type of reads
        - Hybrid: Using a single assembler with both types of read in one run (only hifiasm --ul)
        - Scaffold: Separately assembling ONT and HiFi reads, and then scaffolding one onto the other

    Each sample can only have one strategy, and branching happens here.
    */
    ch_main_branched
            .to_assemble
            .branch { it ->
                single: it.strategy == "single"
                hybrid: it.strategy == "hybrid"
                scaffold: it.strategy == "scaffold"
            }
            .set { ch_main_assemble_branched }

    ch_main_assemble_branched
        .single
        .dump(tag: "Assemble: Branched: Single")
    ch_main_assemble_branched
        .hybrid
        .dump(tag: "Assemble: Branched: Hybrid")
    ch_main_assemble_branched
        .scaffold
        .dump(tag: "Assemble: Branched: scaffold")

    /*
    Inputs for flye assembler:
        - Samples with single strategy, where the assembler is flye
        - Samples from the scaffold strategy where either (or both) assembler is flye
    */
    ch_main_assemble_branched
        .single
        .filter { it -> it.assembler1 == "flye" }
        .mix(
            // Does this actually work correctly? What happens to samples where assembler1 and assembler2 are flye?
            // TODO: May need fixing, think those need to be mixed in individually so they actually are assembled twice
            ch_main_assemble_branched
                .scaffold
                .filter { it -> it.assembler1 == "flye" || it.assembler2 == "flye" }
        )
        .set { ch_main_assemble_flye }

    // Assembly flye branch
    // Extra args per sample are stored in the meta map, so is the estimated / expected genome size
    ch_main_assemble_flye
        .multiMap {
            it ->
            reads: [
                [
                    id: it.meta.id,
                    genome_size: it.genome_size,
                    flye_args: it.flye_args ?: ""
                ],
                // Reads are matched based on assembler
                // Does this actually work correctly? What happens to samples where assembler1 and assembler2 are flye?
                // TODO: May need fixing
                it.assembler1 == "flye" ? it.ontreads : (it.assembler2 == "flye" ? it.hifireads : []),
            ]
            mode: it.assembler1 == "flye" ? "--nano-hq" : "--pacbio-hifi"
        }
        .set { flye_inputs }

    flye_inputs.reads.dump(tag: "Assemble: Flye inputs")

    // Run through flye
    FLYE(flye_inputs.reads, flye_inputs.mode)

    ch_versions = ch_versions.mix(FLYE.out.versions)

    /* Hifiasm: everything that is not hifiasm-ONT
        Single branch with hifiasm as assembler and no ont reads (only hifireads)
        Hybrid assembly
        Scaffold samples where assembler2 (hifi assembler) is hifiasm
    */
    ch_main_assemble_branched
            .single
            .filter { it -> it.assembler1 == "hifiasm" && !it.ontreads }
            .mix(
                ch_main_assemble_branched
                    .hybrid
                    .filter { it -> it.assembler1 == "hifiasm" }
            )
            .mix(ch_main_assemble_branched
                    .scaffold
                    .filter { it -> it.assembler2 == "hifiasm"  }
                    // the samples for scaffolding should not have ONT reads, otherwise hifiasm will run in --ul mode
                    .map { it -> it - it.subMap("ontreads") }
            )
            .set { ch_main_assemble_hifi_hifiasm }

    ch_main_assemble_hifi_hifiasm.dump(tag: "Assemble: hifiasm HIFI inputs")



    HIFIASM(ch_main_assemble_hifi_hifiasm
                .map {
                    it -> [
                        // Put sample-level args into meta map
                        [id: it.meta.id, hifiasm_args: it.hifiasm_args ?: ""],
                        it.hifireads,
                        // for hybrid samples include ONT reads in 3rd slot of first input (see hifiasm module)
                        (it.stragtegy == "hybrid" && it.ontreads) ? it.ontreads : []
                        ]
                    },
            [[], [], []],
            [[], [], []],
            [[], []])

    // hifiasm produces GFA files, convert to fasta & restore meta map with id only
    GFA_2_FA_HIFI( HIFIASM.out.processed_unitigs.map { meta, fasta -> [[id: meta.id], fasta] } )

    ch_versions = ch_versions.mix(HIFIASM.out.versions).mix(GFA_2_FA_HIFI.out.versions)

    /*
    Assemble hifiasm_ont branch:
        Single branch with hifiasm and only ont reads
        Scaffold samples where assembler1 (ont assembler) is hifiasm
    */
    ch_main_assemble_branched
        .single
        .filter { it -> it.assembler1 == "hifiasm" && it.ontreads }
        .mix(ch_main_assemble_branched
                .scaffold
                .filter { it -> it.assembler1 == "hifiasm"  }
        )
        .set { ch_main_assemble_ont_hifiasm }

    ch_main_assemble_ont_hifiasm.dump(tag: "Assemble: hifiasm ONT inputs")

    HIFIASM_ONT(ch_main_assemble_ont_hifiasm.map { it -> [ [id: it.meta.id, hifiasm_args: it.hifiasm_args ?: ""],  it.ontreads, [] ] }, [[], [], []], [[], [], []], [[], []])

    GFA_2_FA_ONT( HIFIASM_ONT.out.processed_unitigs.map { meta, fasta -> [[id: meta.id], fasta] } )

    ch_versions = ch_versions.mix(HIFIASM_ONT.out.versions).mix(GFA_2_FA_ONT.out.versions)


    // Now, the individual assemblies need to be correctly added into the main channel.
    // This should be done per-strategy I think
    // join assembler outputs back to assembler inputs and determine correct placement of the assembly.

    // Flye:
    ch_main_assemble_flye
        // Convert to list for join
        .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        .join( FLYE.out.fasta
                .map { meta, assembly -> [meta: [id: meta.id], flye_assembly: assembly ] }
                .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        )
        // After joining re-create the maps from the stored map
        .map { it -> it.collect { _entry, map -> [ (map.key): map.value ] }.collectEntries() }
        // The flye_assembly has to be placed into the correct slot
        // TODO: Rethink if this works for flye/flye scaffolds?!
        .map { it -> it - it.subMap("flye_assembly") +
                [
                    assembly:  it.strategy == "single" ? it.flye_assembly : null,
                    assembly1: it.assembler1 == "flye" ? it.flye_assembly : null,
                    assembly2: it.assembler2 == "flye" ? it.flye_assembly : null,
                ]
        }
        .set { flye_assemblies }
    flye_assemblies.dump(tag: "Assemble: Flye assemblies")

    // Join hifiasm hifi assemblies back to main channel
    ch_main_assemble_hifi_hifiasm
        // Convert to list for join
        .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        .join( GFA_2_FA_HIFI.out.contigs_fasta
                .map { meta, assembly -> [meta: meta, hifiasm_assembly: assembly ] }
                .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        )
        // After joining re-create the maps from the stored map
        .map { it -> it.collect { _entry, map -> [ (map.key): map.value ] }.collectEntries() }
        // On the proper map, place the hifiasm assembly into the correct position
        .map {
            // remove hifiasm_assembly entry from joined result
            it -> it - it.subMap("hifiasm_assembly") +
            // stick what was in hifiasm_assembly into the correct key
            [
                assembly: (it.strategy == "single" || it.strategy == "hybrid") && it.assembler1 == "hifiasm" ? it.hifiasm_assembly : null,
                // I think below case dose not exist in this channel since it is only hifiasm (assembler2) assemblies?
                assembly1: it.strategy == "scaffold" && it.assembler1 == "hifiasm" ? it.hifiasm_assembly : null,
                assembly2: it.strategy == "scaffold" && it.assembler2 == "hifiasm" ? it.hifiasm_assembly : null
            ]
        }
        .set { hifiasm_hifi_assemblies }

    hifiasm_hifi_assemblies.dump(tag: "Assemble: hifiasm HIFI assemblies")

    // Join hifiasm ONT assemblies back to main channel
    ch_main_assemble_ont_hifiasm
        .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        .join( GFA_2_FA_ONT.out.contigs_fasta
                .map { meta, assembly -> [meta: meta, hifiasm_assembly: assembly ] }
                .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        )
        // After joining re-create the maps from the stored map
        .map { it -> it.collect { _entry, map -> [ (map.key): map.value ] }.collectEntries() }
        .map {
            it -> it -it.subMap("hifiasm_assembly") +
            [
                assembly: (it.strategy == "single" || it.strategy == "hybrid") && it.assembler1 == "hifiasm" ? it.hifiasm_assembly : null,
                // I think below case dose not exist in this channel since it is only ont (assembler1) assemblies?
                assembly1: it.strategy == "scaffold" && it.assembler1 == "hifiasm" ? it.hifiasm_assembly : null,
                assembly2: it.strategy == "scaffold" && it.assembler2 == "hifiasm" ? it.hifiasm_assembly : null
            ]
        }
        .set { hifiasm_ont_assemblies }

    hifiasm_ont_assemblies.dump(tag: "Assemble: hifiasm HIFI assemblies")

    // The single and hybrid channels can be mixed and forwarded.
    // The scaffold channel needs to be joined separately.
    flye_assemblies
        .filter { it -> ["single","hybrid"].contains(it.strategy) }
        .mix(
            hifiasm_hifi_assemblies
                .filter { it -> ["single","hybrid"].contains(it.strategy) }
        )
        .mix(
            hifiasm_ont_assemblies
                .filter { it -> ["single","hybrid"].contains(it.strategy) }
        )
        .set { ch_assemblies_no_scaffold }

    ch_assemblies_no_scaffold.dump(tag: "Assemble: Assemblies without scaffolding")

    // This leaves the scaffold strategy.
    // scaffolds can be: FLYE-HIFIASM, FLYE-FLYE, HIFIASM-HIFIASM HIFIASM-FLYE or

    flye_assemblies
        // Flye-hifiasm
        .filter { it -> it.strategy == "scaffold" && it.assembler1 == "flye" && it.assembler2 == "hifiasm" }
        .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        .join( hifiasm_hifi_assemblies
                .filter{ it -> it.strategy == "scaffold" && it.assembler1 == "flye" && it.assembler2 == "hifiasm" }
                .map { it -> [ meta: it.meta, hifiasm_assembly: it.assembly2 ] }
                .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        )
        .map { it -> it.collect { _entry, map -> [ (map.key): map.value ] }.collectEntries() }
        .map { it -> it - it.subMap("hifiasm_assembly","assembly2") + [assembly2: it.hifiasm_assembly] }
        .set{ scaffold_flye_hifiasm }

    // flye-flye
    flye_assemblies
        .filter { it -> it.strategy == "scaffold" && it.assembler1 == "flye" && it.assembler2 == "flye" }
        .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        .join( flye_assemblies
                .filter{ it -> it.strategy == "scaffold" && it.assembler1 == "flye" && it.assembler2 == "flye" }
                .map { it -> [ meta: it.meta, flye_assembly: it.assembly2 ] }
                .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        )
        .map { it -> it.collect { _entry, map -> [ (map.key): map.value ] }.collectEntries() }
        .map { it -> it - it.subMap("flye_assembly", "assembly2") + [assembly2: it.flye_assembly] }
        .set{ scaffold_flye_flye }

    // hifiasm_flye
    hifiasm_ont_assemblies
        .filter { it -> it.strategy == "scaffold" && it.assembler1 == "hifiasm" && it.assembler2 == "flye" }
        .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        .join( flye_assemblies
                .filter{ it -> it.strategy == "scaffold" && it.assembler1 == "hifiasm" && it.assembler2 == "flye" }
                .map { it -> [ meta: it.meta, flye_assembly: it.assembly2 ] }
                .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        )
        .map { it -> it.collect { _entry, map -> [ (map.key): map.value ] }.collectEntries() }
        .map { it -> it - it.subMap("flye_assembly", "assembly2") + [assembly2: it.flye_assembly] }
        .set{ scaffold_hifiasm_flye }

    // hifiasm_hifiasm
    hifiasm_ont_assemblies
        .filter { it -> it.strategy == "scaffold" && it.assembler1 == "hifiasm" && it.assembler2 == "hifiasm" }
        .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        .join( hifiasm_hifi_assemblies
                .filter{ it -> it.strategy == "scaffold" && it.assembler1 == "hifiasm" && it.assembler2 == "hifiasm" }
                .map { it -> [ meta: it.meta, hifiasm_assembly: it.assembly2 ] }
                .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        )
        .map { it -> it.collect { _entry, map -> [ (map.key): map.value ] }.collectEntries() }
        .map { it -> it - it.subMap("hifiasm_assembly","assembly2") + [assembly2: it.hifiasm_assembly] }
        .set{ scaffold_hifiasm_hifiasm }

    // branch to scaffold those assemblies that need it

    scaffold_flye_hifiasm
        .mix(scaffold_flye_flye)
        .mix(scaffold_hifiasm_flye)
        .mix(scaffold_hifiasm_hifiasm)
        .set { ch_to_scaffold }

    ch_to_scaffold.dump(tag: "Assemble: Assemblies with scaffolding - inputs")

    // For scaffolding, depeding on which strategy we used, the correct assembly needs to go into either target or query:
    // assembly1 is always ONT, assembly2 is always HiFi

    ch_to_scaffold
        .multiMap {
            it ->
            target: [
                it.meta,
                it.assembly_scaffolding_order == "ont_on_hifi" ? (it.assembly1) : (it.assembly2)
                ]
            query: [
                it.meta,
                it.assembly_scaffolding_order == "ont_on_hifi" ? (it.assembly2) : (it.assembly1)
                ]
        }
        .set { ragtag_in }

    // Scaffold with PATCH
    RAGTAG_PATCH(ragtag_in.target, ragtag_in.query, [[], []], [[], []] )

    // Update inputs
    ch_to_scaffold
        .map { it -> it - it.subMap("assembly") }
        .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        .join(
            RAGTAG_PATCH.out.patch_fasta
                .map { it -> [meta: it[0], assembly: it[1]] }
                .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        )
        .map { it -> it.collect { _entry, map -> [ (map.key): map.value ] }.collectEntries() }
        .set { ch_assemblies_scaffold }

    ch_assemblies_scaffold.dump(tag: "Assemble: Assemblies with scaffolding - outputs")

    // Mix everything assembled back togehter
    ch_assemblies_no_scaffold
        .mix(ch_assemblies_scaffold)
        .set { ch_main_assembled }

    ch_main_assembled.dump(tag: "Assemble: Assembled")

    ch_versions = ch_versions.mix(RAGTAG_PATCH.out.versions)

    // Mix with whatever was not destined for assembly
    ch_main_branched
        .no_assemble
        .mix( ch_main_assembled )
        .set { ch_main_to_mapping }

    ch_main_to_mapping.dump(tag: "Assemble: TO MAPPING")


    ch_main_to_mapping
        .branch {
            it ->
            quast: it.quast
            no_quast: !it.quast
        }
        // Note that this channel is set here but the quast branch is further used
        .set { ch_main_quast_branch }

    // If QUAST should run, and we need an alignment to reference, this is created here
    ch_main_quast_branch
        .quast
        .branch {
            it ->
                use_ref: it.use_ref
                no_use_ref: !it.use_ref
        }
        .set {
            ch_quast_branched
        }
    // It is actually only created if no bam file is provided
    ch_quast_branched
        .use_ref
        .branch { it ->
            to_map: !it.ref_map_bam
            dont_map: it.ref_map_bam
        }
        .set { ch_ref_mapping_branched }

    // Use the QC reads and map them to ref
    ch_ref_mapping_branched
        .to_map
        .map {
            it ->
            [ [id: it.meta.id, qc_reads: it.qc_reads], it.qc_reads_path, it.ref_fasta ]
        }
        .set { map_to_ref_in }

    MAP_TO_REF(map_to_ref_in) // returns meta: [ id: ]

    // Add the ref mapping to the large main channel
    ch_ref_mapping_branched
        .to_map
        .map { it -> it - it.subMap("ref_map_bam") }
        .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        .join(
            MAP_TO_REF.out.ch_aln_to_ref_bam
            // Note that this is a normal list channel and needs to become a map before conversion back to list and joining
            // Otherwise the map cannot be regenerated later
                .map { it -> [meta: it[0], ref_map_bam: it[1]] }
                .map { it -> it.collect { entry -> [ entry.value, entry ] } }
        )
        .map { it -> it.collect { _entry, map -> [ (map.key): map.value ] }.collectEntries() }
        .mix(ch_ref_mapping_branched.dont_map)
        .mix(ch_quast_branched.no_use_ref)
        // above recreates ch_main_quast_branch.quast
        .mix(ch_main_quast_branch.no_quast)
        .set { ch_main_to_qc }


    //QC on initial assembly

    // scaffolds to QC need to be defined here, this is what is in the assembly slot
    ch_main_to_qc
        .map { it -> [it.meta, it.assembly] }
        .set { scaffolds }

    QC(ch_main_to_qc, scaffolds, meryl_kmers)

    ch_versions = ch_versions.mix(QC.out.versions)

    // If annotation liftover on the initial assembly is desired, it happens here.
    ch_main_to_qc
        .filter {
            it -> it.lift_annotations
        }
        .map { it ->
            [
                it.meta,
                it.assembly,
                it.ref_fasta,
                it.ref_gff
            ]
        }
        .set { liftoff_in }

    RUN_LIFTOFF(liftoff_in)
    ch_versions = ch_versions.mix(RUN_LIFTOFF.out.versions)

    emit:
    ch_main                     = ch_main_to_qc
    assembly_quast_reports      = QC.out.quast_out
    assembly_busco_reports      = QC.out.busco_out
    assembly_merqury_reports    = QC.out.merqury_report_files
    versions                    = ch_versions
}
