include { DOWNLOAD_REF } from './alignment/download-req.nf'
include { BWA_MEM } from './alignment/bwa_alignment.nf'
include { HLA_LA } from './HLA-LA/main.nf'
include { HLA_HD } from './HLA-HD/main.nf'
include { XHLA } from './xHLA/main.nf'

workflow MAPPING {
    main:
        if (params.ref && params.refWoAlt) {
            references = Channel.fromPath([
                    params.ref + "*",
                    params.refWoAlt + "*"]
                ).buffer(size: 6)
        } else {
            DOWNLOAD_REF()
            references = DOWNLOAD_REF.out.ref.concat(DOWNLOAD_REF.out.refWoAlt)
        }
        fa_files = Channel.fromFilePairs('1-Input/*R{1,2}*.fastq')
        BWA_MEM(fa_files, references)
    emit:
        ref = BWA_MEM.out.first()
        refWoAlt = BWA_MEM.out.last()
}

workflow HLATYPING_REF {
    take: ref
    main:
        HLA_LA(ref)
        HLA_HD(ref)
}

workflow HLATYPING_REF_WO_ALT {
    take: ref
    main:
        XHLA(ref)
}

workflow {
    MAPPING()
    HLATYPING_REF(MAPPING.out.ref)
    HLATYPING_REF_WO_ALT(MAPPING.out.refWoAlt)
}