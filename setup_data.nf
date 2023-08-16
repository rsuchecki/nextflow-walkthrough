#!/usr/bin/env nextflow

nextflow.enable.dsl=2 

//Take accessions defined in nextflow.config.
//Use --take N to process first N accessions or --take all to process all
ACCESSIONS_CHANNEL = Channel.from(params.accessions).take( params.download == 'all' ? -1 : params.download )

region = "${params.chr}_${params.start}-${params.end}"

process GET_ADAPTERS { //alternative: referencesChannel = Channel.fromPath(params.reference)
  publishDir "data/misc", mode: 'copy'

  input: 
    val(adapters_url)

  output:
    path('*')

  script:
  """
  wget ${adapters_url}
  """
}

process GET_REFERENCE { //alternative: referencesChannel = Channel.fromPath(params.reference)
  tag { "${region}"}
  publishDir "data/references", mode: 'copy'

  input:
    val(reference_url)

  output:
    path('*')

  script:
  """
  wget ${reference_url}
  """
}

process GET_READS {
  tag { "${accession} @ ${region}"}
  publishDir "data/raw_reads", mode: 'copy'

  input:
    val(accession)
    //e.g. ACBarrie

  output:
    tuple val(accession), path("${accession}_R?.fastq.gz")
    //e.g. ACBarrie, [ACBarrie_R1.fastq.gz, ACBarrie_R2.fastq.gz]

  script:
  URL_BASE = [params.reads_base_url, region, accession].join('/')
  """
  wget ${URL_BASE}_${params.r1_suffix} ${URL_BASE}_${params.r2_suffix} \
  && zcat ${accession}_R1.fastq.gz | head | awk 'END{exit(NR<4)}' \
  && zcat ${accession}_R2.fastq.gz | head | awk 'END{exit(NR<4)}'
  """
}

workflow {
  //Reference fasta.gz specified in nextflow.config, override with --reference url_to_ref.fasta.gz
  GET_REFERENCE(params.reference)
  
  //Adapters FASTA file URL specified in nextflow.config, override with --adapters url_to_adapters.fasta
  GET_ADAPTERS(params.adapters)
  
  //
  GET_READS(ACCESSIONS_CHANNEL)
}