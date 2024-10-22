include { BUSCO } from '../../../../modules/local/busco/main'

workflow RUN_BUSCO {
  take: 
    assembly

  main:
    if(params.busco) BUSCO(assembly, params.busco_lineage, params.busoc_db)
}