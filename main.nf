#!/usr/bin/env nextflow

// DSL1 is used by default, switching to DSL2
nextflow.enable.dsl=2

/*
  process 1 file unless --n used at run time, e.g. --n 32 
  to process all FASTQ files (16 pairs)
*/
params.n = 1 

Channel.fromPath("data/raw_reads/*.fastq.gz")
  .take( params.n )
  .set { ReadsForQcChannel }

process FASTQC {  
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

Channel.fromPath('data/references/reference.fasta.gz')
  .set { ReferencesChannel }

process BWA_INDEX {
  input:
    path(ref)

  output:
    tuple val("${ref}"), path("*") 

  script:
  """
  bwa index -a bwtsw ${ref}
  """
}

Channel.fromFilePairs("data/raw_reads/*_R{1,2}.fastq.gz")
  .take( params.n )
  .set{ ReadPairsForTrimmingChannel }

Channel.fromPath('data/misc/TruSeq3-PE.fa')
.set{ AdaptersChannel }

process TRIM_PE {
  tag { "$sample" }
  input:
    tuple val(sample), path(reads), path(adapters) 

  output:
    tuple val(sample), path('*.paired.fastq.gz') 

  script:
  """
  trimmomatic PE \
  ${reads} \
  R1.paired.fastq.gz \
  R1.unpaired.fastq.gz \
  R2.paired.fastq.gz \
  R2.unpaired.fastq.gz \
  ILLUMINACLIP:${adapters}:2:30:10:3:true \
  LEADING:2 \
  TRAILING:2 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36 \
  -Xms256m \
  -Xmx256m
  """
}

/*
 Chaining everything toogether
*/
workflow {
  //QC - could be separated as a sub-workflow
  FASTQC( ReadsForQcChannel )
  MULTIQC( FASTQC.out.collect() )

  //Workflow proper
  BWA_INDEX(  ReferencesChannel )
  TRIM_PE ( ReadPairsForTrimmingChannel.combine( AdaptersChannel ) )
}