process BWA_MEM {
    cpus 8
    memory '20G'
    container 'docker://humanlongevity/hla'
    publishDir "2-Alignment", mode: 'link'
    input:
        tuple val(sample), path(reads)
        path bwa_index
    output:
        tuple val(sample), path("${sample}.bam*")
    script:
        def idxbase = bwa_index[0].baseName
        def (f1, f2) = reads
        """
        bwa mem -t ${task.cpus} ${idxbase} ${f1} ${f2} |\
        samtools sort -@${task.cpus} -o ${sample}.bam -
        samtools index ${sample}.bam
        """
}