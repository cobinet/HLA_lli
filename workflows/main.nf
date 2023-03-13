include { DOWNLOAD_REF } from './alignment/download-req.nf'
include { BWA_MEM } from './alignment/bwa_alignment.nf'
include { HLA_LA } from './HLA-LA/main.nf'
include { HLA_HD } from './HLA-HD/main.nf'
include { XHLA } from './xHLA/main.nf'
include { OPTITYPE } from './Optitype/main.nf'
include { EXTRACT_MHC } from './utils/main.nf'
include { BAM_TO_FQ } from './utils/main.nf'

workflow MAPPING {
    main:
        // Get reads
        fa_files = Channel.fromFilePairs('1-Input/*R{1,2}*.fastq')
        // Get reference genome
        if (params.ref && params.refWoAlt) {
            references = Channel.fromPath([
                    params.ref + "*",
                    params.refWoAlt + "*"]
                ).buffer(size: 6)
        } else {
            DOWNLOAD_REF()
            references = DOWNLOAD_REF.out.ref.concat(DOWNLOAD_REF.out.refWoAlt)
        }
        // Map reads to reference
        BWA_MEM(fa_files, references)
    emit:
        ref = BWA_MEM.out.first()
        refWoAlt = BWA_MEM.out.last()
}

workflow HLATYPING_REF {
    take: bam
    main:
        HLA_LA(bam)

        filtered_reads = EXTRACT_MHC(bam) | BAM_TO_FQ
        HLA_HD(filtered_reads)
        OPTITYPE(filtered_reads)
}

workflow HLATYPING_REF_WO_ALT {
    take: bam
    main:
        XHLA(bam)
}

workflow {
    MAPPING()
    HLATYPING_REF(MAPPING.out.ref)
    HLATYPING_REF_WO_ALT(MAPPING.out.refWoAlt)
}