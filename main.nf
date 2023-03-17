#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HLATYPING } from './workflows/main.nf'

workflow {
    fq_files = Channel.fromFilePairs('1-Input/*R{1,2}*.fastq')
    HLATYPING(fq_files)
}
