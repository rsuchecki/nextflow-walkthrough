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
