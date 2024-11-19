include { PILON } from '../../../../../modules/nf-core/pilon/main'
workflow RUN_PILON {
    take:
      assembly_in
      aln_to_assembly_bam_bai

    main:
      assembly_in
        .join(aln_to_assembly_bam_bai)
        .set { pilon_in }

      PILON(pilon_in.map {it-> [it[0],it[1],it[2]]},
            pilon_in.map {it-> [it[0],it[3],it[4]]}, 
            "bam")
    
    emit:
      PILON.out.improved_assembly
}