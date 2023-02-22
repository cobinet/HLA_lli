include { DOWNLOAD_REF } from './alignment/download-req.nf'
include { BWA_MEM } from './alignment/bwa_alignment.nf'


workflow HLATYPING {
    emit:
        BWA_MEM.out
    main:
        DOWNLOAD_REF()
        fa_files = Channel.of.fromFilePairs('1-Input/*R{1,2}*.fastq')
        BWA_MEM(fa_files , DOWNLOAD_REF.out)
}

workflow {
    HLATYPING()
}