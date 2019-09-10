Channel.fromPath("data/raw_reads/*.fastq.gz")
  .take(1)
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
