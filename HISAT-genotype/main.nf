
process hisat_genotype {
    cpus 8
    memory '40G'
    container 'docker://nmendozam/hisat-genotype:3.50.0'
    publishDir "2-results", mode: 'copy'
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("hisatgenotype_out/*.report")
    script:
        def (f1, f2) = reads
        """
        hisatgenotype --base hla -1 ${f1} -2 ${f2}
        """
}

workflow {
    reads = Channel
    .fromFilePairs('1-Input/*R{1,2}*.fastq')
    // .fromFilePairs("1-Input/J27462-S1-L3_S1_L001_R{1,2}_001.fastq")

    hisat_genotype(reads)
}
