process download_reference {
    container 'frolvlad/alpine-bash'
    output:
        path "hg38.fa", emit: "reference"
    script:
        """
        wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
        gunzip hg38.fa.gz
        """
}

process remove_alternatives {
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
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
    container 'biocontainers/bwa:v0.7.17_cv1'
    publishDir "reference_genome"
    input:
        file(fasta)
    output:
        path "${fasta}.{,amb,ann,bwt,pac,sa}"
    script:
        """
        bwa index ${fasta}
        """
}

workflow DOWNLOAD_REF {
    main:
        ref = download_reference()
        noalt = remove_alternatives(ref)

        index_fasta(ref.concat(noalt))

    emit:
        ref = index_fasta.out[0]
        refWoAlt = index_fasta.out[1]
}

workflow {
    DOWNLOAD_REF()
}
