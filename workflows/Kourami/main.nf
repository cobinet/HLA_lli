process TYPING {
    container "aniroula/kourami:latest"
    input:
        tuple val(sample), path(files)
    output:
        path "out"

    script:
        script_path = "/usr/local/kourami/scripts"
        target_path =  "/usr/local/kourami/target"
        db = "/usr/local/kourami/db/"
        def (bam, bai) = files
        """
        $script_path/alignAndExtract_hs38DH.sh -d $db $sample $bam
        java -jar $target_path/Kourami.jar -d $db ${sample}_on_KouramiPanel.bam
        """
}

workflow KOURAMI {
    take: align_files
    main:
        TYPING(align_files)
    emit:
        TYPING.out
}

workflow {
    align_files = Channel.fromFilePairs('alignment/*.{bam,bam.bai}')
    KOURAMI(align_files)
}