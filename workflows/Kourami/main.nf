process MAKE_DB {
    container "aniroula/kourami:latest"
    input:
        val ulr "https://github.com/ANHIG/IMGTHLA/archive/refs/tags/v3.50.0-alpha.tar.gz"
    output:
        path "db"
    script:
        """
        /usr/local/kourami/scripts/download_panel.sh

        wget -O- $url | tar xzf - -C db
        /usr/local/kourami/scripts/formatIMGT.sh
        bwa index db/All_FINAL_with_Decoy.fa.gz
        """
}
process KOURAMI {
    container "aniroula/kourami:latest"
    input:
        tuple val(sample), path(bam)
    output:

    script:
        script_path = "/usr/local/kourami/scripts"
        target_path =  "/usr/local/kourami/target"
        db = "/usr/local/kourami/db/"
        """
        $script_path/alignAndExtract_hs38DH.sh -d $db $sample $bam
        java -jar $target_path/Kourami.jar -d $db ${sample}_on_KouramiPanel.bam
        """
}