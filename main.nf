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
    tuple  val(sample), path(reads), path(adapters) 

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


process BWA_ALIGN {
  tag { "$sample" }
  publishDir 'results/aligned', mode: 'copy'

  input:
    tuple val(sample), path(reads), val(prefix), path(index) 

  output:
    path '*.bam'

  script:
  """
  bwa mem -t ${task.cpus} ${prefix} ${reads} \
  | samtools view -b > ${sample}.bam
  """
}

process MERGE_BAMS {
  publishDir 'results/merged', mode: 'copy'
  cpus 2

  input:
    path(BAMs) // from alignedReadsChannel.collect()

  output:
    path('*.bam') //Input BAMs will be omitted, could also be explicit:  path("${params.n}_samples_megred.bam")

  script:
  """
  samtools merge --threads ${task.cpus} ${params.n}_samples_megred.bam *.bam
  """
}

/*
 Chaining everything toogether
*/
workflow {
  /*
  For demonstration purposes we are comparing 3 different 
  syntax styles of workflow composition
  These can also be mixed. 
  The 'pipes' work best when each process has only one input (which may be a tuple of multiple elements)
  */
  if(params.pipes) {                                              //Syntax style 1
    //QC - could be separated as a sub-workflow
    ReadsForQcChannel | FASTQC | collect | MULTIQC
    
    //Workflow proper    
    ReadPairsForTrimmingChannel \
    | combine( AdaptersChannel ) \
    | TRIM_PE \
    | combine( ReferencesChannel | BWA_INDEX ) \
    | BWA_ALIGN \
    | collect \
    | MERGE_BAMS

  } else if(params.nested) {                                       //Syntax style 2
    //QC - could be separated as a sub-workflow
    MULTIQC ( FASTQC( ReadsForQcChannel ).collect() )

    //Workflow proper
    MERGE_BAMS (
      BWA_ALIGN (      
        TRIM_PE ( 
          ReadPairsForTrimmingChannel
          .combine( AdaptersChannel ) 
        )
        .combine( BWA_INDEX( ReferencesChannel ) )
      )
      .collect()
    )
  } else {                                                         //Syntax style 3
    //QC - could be separated as a sub-workflow
    FASTQC( ReadsForQcChannel )
    MULTIQC( FASTQC.out.collect() )

    //Workflow proper
    TRIM_PE ( ReadPairsForTrimmingChannel.combine( AdaptersChannel ) )
    BWA_INDEX(  ReferencesChannel )
    BWA_ALIGN ( TRIM_PE.out.combine( BWA_INDEX.out ) )
    MERGE_BAMS ( BWA_ALIGN.out.collect() )    
  }
}