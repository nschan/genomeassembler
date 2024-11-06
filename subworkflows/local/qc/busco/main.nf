include { BUSCO } from '../../../../modules/local/busco/main'

workflow RUN_BUSCO {
  take: 
    assembly

  main:
  Channel.empty().set { busco_batch_summary }
  Channel.empty().set { busco_short_summary_txt }
  Channel.empty().set { busco_short_summary_json }

  if(params.busco) {
      BUSCO(assembly, params.busco_lineage, params.busoc_db ? file( params.busoc_db, checkIfExists: true ) : [])
      BUSCO
        .out
        .batch_summary
        .set { busco_batch_summary }
      BUSCO
        .out
        .short_summaries_txt
        .set { busco_short_summary_txt }
      BUSCO
        .out
        .short_summaries_json
        .set { busco_short_summary_json }
  }
  
  emit:
    busco_batch_summary
    busco_short_summary_json
    busco_short_summary_txt
}