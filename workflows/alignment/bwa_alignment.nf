process BWA_MEM {
    cpus 8
    memory '20G'
    container 'docker://humanlongevity/hla'
    publishDir "2-Alignment"
    input:
        tuple val(sample), path(reads)
        each path(bwa_index)
    output:
        tuple val(sample), path("${sample}.*.bam*")
    script:
        def idxbase = bwa_index[0].baseName
        def refname = bwa_index[0].simpleName
        def (f1, f2) = reads
        """
        bwa mem -t ${task.cpus} ${idxbase} ${f1} ${f2} |\
        samtools sort -@${task.cpus} -o ${sample}.${refname}.bam -
        samtools index ${sample}.${refname}.bam
        """
}