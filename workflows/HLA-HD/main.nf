params.hlahd_dir = "~/.local/share/hlahd.1.6.1/"
params.reference_ver = "GRCh38"

process TYPING {
    memory '10G'
    cpus '20'
    container "nmendozam/hla-hd"
    publishDir "2-HLA-HD", mode: 'move'
    input:
        tuple val(sample), path(reads)
    output:
        path "${sample}/result/*"
    script:
        def (f1, f2) = reads
        """
        hlahd.sh -t 20 -m 100 -c 0.95 \
        -f ${params.hlahd_dir}freq_data/ \
        ${f1} ${f2} \
        ${params.hlahd_dir}HLA_gene.split.txt \
        ${params.hlahd_dir}dictionary/ \
        ${sample} \
        .
        """

}

process uncompress_fastq {
    publishDir "2-Uncompresed", mode: 'move'
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_{1,2}.fastq")
    script:
        def (f1, f2) = reads
        """
        gzip -cd ${f1} > ${sample}_1.fastq
        gzip -cd ${f2} > ${sample}_2.fastq
        """
}

process EXTRACT_MHC {
    container "biocontainers/samtools:v1.9-4-deb_cv1"
    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path("${sample}.mhc.bam")
    script:
        if (params.reference_ver == "GRCh38") {
            region = "chr6:28,510,120-33,480,577"
        }
        else if (params.reference_ver == "GRCh37"){
            region = "chr6:28,477,797-33,448,354"
        }
        else {
            error "Expected reference version either GRCh38 or GRCh37"
        }
        """
        #Extract MHC region
        samtools view -h -b ${bam} ${region} > ${sample}.mhc.bam
        
        #Extract unmap reads
        samtools view -b -f 4 ${bam} > ${sample}.unmap.bam
        
        #Merge bam files
        samtools merge -o ${sample}.merge.bam ${sample}.unmap.bam ${sample}.mhc.bam
        """
}

process BAM_TO_FQ {
    container "broadinstitute/picard"
    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path("${sample}.hla.{1,2}.fastq")
    script:
        """
        #Create fastq
        java -jar picard.jar SamToFastq I=${bam} F=${sample}.hlatmp.1.fastq F2=${sample}.hlatmp.2.fastq
        
        #Change fastq ID
        cat ${sample}.hlatmp.1.fastq |awk '{if(NR%4 == 1){O=\$0;gsub("/1","1",O);print O}else{print \$0}}' > ${sample}.hla.1.fastq
        cat ${sample}.hlatmp.2.fastq |awk '{if(NR%4 == 1){O=\$0;gsub("/2","2",O);print O}else{print \$0}}' > ${sample}.hla.2.fastq
        """
}

workflow HLA_HD {
    take:
        bam
    main:
        EXTRACT_MHC(bam) | BAM_TO_FQ | TYPING
    emit:
        TYPING.out
}

workflow {
    reads = Channel
    .fromFilePairs('1-Input/*R{1,2}*.fastq')

    hla_hd(reads)
}
