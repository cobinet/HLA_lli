process download_reference {
    container 'docker://alpine'
    output:
        path "hg38.fa", emit: "reference"
    script:
        """
        wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
        gunzip hg38.fa.gz
        """
}

process remove_alternatives {
    container 'docker://alpine'
    input:
        file reference
    script:
        def out_file = "${reference.baseName}-noalt.fa"
        """
        samtools faidx -o ${out_file} ${reference}\
            `grep "^>" ${reference} | grep -v "alt" | cut -c 2-`
        """
    output:
        path "*-noalt.fa"
}

process index_fasta {
    input:
        file(fasta)
    output:
        path "${fasta}.{,amb,ann,bwt,pac,sa}"
    script:
        """
        bwa index ${fasta}
        """
}

workflow {

    remove_alternatives(download_reference()) | index_fasta
}