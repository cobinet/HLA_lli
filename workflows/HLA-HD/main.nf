process HLA_HD {
    tag "${sample}"
    memory '50G'
    cpus '20'
    container "nmendozam/hla-hd:latest"
    publishDir "2-HLA-HD"
    input:
        tuple val(sample), path(reads)
    output:
        path "${sample}/result/*"
    script:
        def (f1, f2) = reads
        hlahd_dir = "/usr/local/hlahd.1.7.0/"
        min_len = 100
        """
        hlahd.sh -t $task.cpus -m $min_len -c 0.95 \
        -f ${hlahd_dir}freq_data/ \
        ${f1} ${f2} \
        ${hlahd_dir}HLA_gene.split.txt \
        ${hlahd_dir}dictionary/ \
        ${sample} \
        .
        """

}

workflow {
    reads = Channel
    .fromFilePairs('1-Input/*R{1,2}*.fastq')

    HLA_HD(reads)
}
