params.hlahd_dir = "~/.local/share/hlahd.1.6.1/"

process hla_hd {
    memory '10G'
    cpus '20'
    // module = ['bowtie2/2.3.5.1']
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

workflow {
    reads = Channel
    .fromFilePairs('1-Input/*R{1,2}*.fastq')

    hla_hd(reads)
}
