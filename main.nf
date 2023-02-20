#!/usr/bin/env nextflow
nextflow.enable.dsl = 2



include { HLATYPING } from './workflows/main'


workflow {
    HLATYPING ()
}
