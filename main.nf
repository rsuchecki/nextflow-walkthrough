params.n = 1 //process 1 file unless --n used at run time, e.g. -n 32

Channel.fromPath("data/raw_reads/*.fastq.gz")
  .take( params.n )
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


Channel.fromFilePairs("data/raw_reads/*_R{1,2}.fastq.gz")
  .take( params.n )
  .set{ readPairsForTrimmingChannel }

Channel.fromPath('data/misc/trimmomatic_adapters/TruSeq3-PE.fa')
.set{ adaptersChannel }

process trimmomatic_pe {
  input:
    set file(adapters), val(accession), file(reads) from adaptersChannel.combine(readPairsForTrimmingChannel)

  output:
    set val(accession), file('*.paired.fastq.gz') into trimmedReadsChannel

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