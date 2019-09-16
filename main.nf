params.n = 1 //process 1 file unless --n used at run time, e.g. -n 32

Channel.fromPath("data/raw_reads/*.fastq.gz")
  .take( params.n )
  .set { readsForQcChannel }

process fastqc {
  input:
    file reads from readsForQcChannel

  """
  fastqc \
    --threads ${task.cpus} \
    ${reads}
  """
}


Channel.fromPath('data/references/reference.fasta.gz')
  .set { referencesChannel }

process bwa_index {
  input:
    file(ref) from referencesChannel

  output:
    set val("${ref}"), file("*") into indexChannel

  script:
  """
  bwa index -a bwtsw ${ref}
  """
}