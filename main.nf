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
  label 'bwa'

  input:
    path(ref)

  output:
    path '*'

  script:
  """
  bwa index -a bwtsw ${ref}
  """
}

Channel.fromFilePairs("data/raw_reads/*_R{1,2}.fastq.gz")
  .take( params.n )
  .set{ ReadPairsForTrimmingChannel }

// If we wanted to treat adapters file as other inputs and pull it down a channel...
// Channel.fromPath('data/misc/TruSeq3-PE.fa').set{ AdaptersChannel }

process TRIM_PE {
  tag { "$sample" }
  input:
    tuple val(sample), path(reads)
    path(adapters) 

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
  MINLEN:36 
  """
}

process BWA_ALIGN {
  label 'align'

  tag { "$sample" }
  publishDir 'results/aligned', mode: 'copy'

  input:
    tuple val(sample), path(reads), path(index) 

  output:
    path '*.bam'

  script:
  """
  bwa mem -t ${task.cpus} ${index[0].baseName} ${reads} \
  | samtools view -b > ${sample}.bam
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
  TRIM_PE ( ReadPairsForTrimmingChannel, file(params.adapters_local) )
  BWA_INDEX (  ReferencesChannel )
  BWA_ALIGN ( TRIM_PE.out.combine(BWA_INDEX.out.toList()) ) 
}