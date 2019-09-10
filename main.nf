params.subset = 1

Channel.fromPath("data/raw_reads/*.fastq.gz")
  .take( params.subset )
  .set { readsForQcChannel }

process fastqc{
  input:
    file reads from readsForQcChannel

  """
  fastqc \
    --threads ${task.cpus} \
    ${reads}
  """
}


referencesChannel = Channel.fromPath('data/references/reference.fasta.gz')

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