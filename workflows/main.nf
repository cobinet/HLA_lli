include { DOWNLOAD_REF } from './alignment/download-req.nf'
include { BWA_MEM } from './alignment/bwa_alignment.nf'

index_pattern = ".{,amb,ann,bwt,pac,sa}"

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
        fa_files = Channel.fromFilePairs('1-Input/*{1,2}_fished.fastq.gz')
        BWA_MEM(fa_files , references)
    emit:
        ref = BWA_MEM.out[0]
        refWoAlt = BWA_MEM.out[1]
}

workflow HLATYPING_REF {
    take:
        path ref
    main:

    emit:

}

workflow HLATYPING_REF_WO_ALT {
    take:
        path ref
    main:

    emit:

}

workflow {
    MAPPING().ref.view()
}