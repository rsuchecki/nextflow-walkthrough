#!/usr/bin/env nextflow

// DSL1 is used by default, switching to DSL2
nextflow.enable.dsl=2

/*
  process 1 file unless --n used at run time, e.g. --n 32 
  to process all FASTQ files (16 pairs)
*/
params.n = 1 
