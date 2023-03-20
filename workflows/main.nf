include { DOWNLOAD_REF } from './alignment/download-req.nf'
include { BWA_MEM } from './alignment/bwa_alignment.nf'
include { HLA_LA } from './HLA-LA/main.nf'
include { HLA_HD } from './HLA-HD/main.nf'
include { XHLA } from './xHLA/main.nf'
include { OPTITYPE } from './Optitype/main.nf'
include { EXTRACT_MHC } from './utils/main.nf'
include { BAM_TO_FQ } from './utils/main.nf'
include { PHLAT } from './PHLAT/main.nf'

workflow GET_REF {
    main:
        // Get reference genome
        if (params.ref && params.refWoAlt) {
            references = Channel.fromPath([
                    params.ref + "*",
                    params.refWoAlt + "*"]
                ).buffer(size: 4)

            if (!references.count()) {
                exit 1, "Wrong reference path"
            }
        } else {
            DOWNLOAD_REF()
            references = DOWNLOAD_REF.out.ref.concat(DOWNLOAD_REF.out.refWoAlt)
        }

    emit:
        references

}

workflow MAPPING {
    take: fq_files
    main:
        // Map reads to reference
        BWA_MEM(fq_files, GET_REF())
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
        PHLAT(filtered_reads)
}

workflow HLATYPING_REF_WO_ALT {
    take: bam
    main:
        XHLA(bam)
}

workflow HLATYPING {
    take: fq_files
    main:
        MAPPING(fq_files)
        HLATYPING_REF(MAPPING.out.ref)
        HLATYPING_REF_WO_ALT(MAPPING.out.refWoAlt)
}

workflow {
    // Get reads
    fq_files = Channel.fromFilePairs('1-Input/*R{1,2}*.fastq')
    HLATYPING(fq_files)
}