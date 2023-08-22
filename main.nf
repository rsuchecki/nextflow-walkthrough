#!/usr/bin/env nextflow

//DSL2 is the default now anyway
nextflow.enable.dsl=2

/*
  process 1 file unless --n used at run time, e.g. --n 16 
  to process all FASTQ files (16 pairs)
*/


Channel.fromPath("data/raw_reads/*.fastq.gz")
  .take( params.n * 2 )
  .set { ReadsForQcChannel }

process FASTQC {  
  module = 'fastqc/0.11.9'
  input:
    path(reads)

  output:
    path('*')

  """
  fastqc \
    --threads ${task.cpus} \
    ${reads}
  """
}

process MULTIQC {
  module = 'multiqc/1.15'
  publishDir 'results/multiqc', mode: 'copy'

  input:
    path('*')

  output:
    path('*')    

  script:
  """
  multiqc .
  """
}

/*
 Chaining everything toogether
*/
workflow {
  //QC - could be separated as a sub-workflow
  FASTQC( ReadsForQcChannel )
  MULTIQC( FASTQC.out.collect() )
}