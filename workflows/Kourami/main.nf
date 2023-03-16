process MAKE_DB {
    container "aniroula/kourami:latest"
    containerOptions "-B db:/usr/local/kourami/db -B custom_db:/usr/local/kourami/custom_db"
    beforeScript "mkdir db custom_db"
    output:
        path "db"
    script:
        def url = "https://github.com/ANHIG/IMGTHLA/archive/refs/tags/v3.50.0-alpha.tar.gz"
        """
        /usr/local/kourami/scripts/download_panel.sh

        wget -O- $url | tar xzf - -C db
        mv db/IMGTHLA-3.50.0-alpha/* db/
        mv db/wmda/* db/
        mv db/alignments/* db/
        mv db/msf/* db/
        /usr/local/kourami/scripts/formatIMGT.sh -i db/ -v 3.50.0 
        bwa index db/All_FINAL_with_Decoy.fa.gz
        """
}

process TYPING {
    container "aniroula/kourami:latest"
    containerOptions "-B ${db}:/usr/local/kourami/${db}"
    input:
        tuple val(sample), path(files)
        path db
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
        TYPING(align_files, MAKE_DB())
    emit:
        TYPING.out
}

workflow {
    align_files = Channel.fromFilePairs('alignment/*.{bam,bam.bai}')
    KOURAMI(align_files)
}