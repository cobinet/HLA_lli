include { DOWNLOAD_REF } from './alignment/download-req.nf'
include { BWA_MEM } from './alignment/bwa_alignment.nf'


workflow HLATYPING {
    main:
        DOWNLOAD_REF()
        fa_files = Channel.fromFilePairs('1-Input/*R{1,2}*.fastq')
        BWA_MEM(fa_files , DOWNLOAD_REF.out.ref)
    emit:
        BWA_MEM.out
}

workflow {
    HLATYPING()
}
