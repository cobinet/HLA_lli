process BWAmem {
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

process xhla {
    cpus 20
    memory '500G'
    container 'docker://humanlongevity/hla'
    publishDir "3-xHLA", mode: 'link'
    input:
        tuple val(sample), file(bam)
    output:
        path "${sample}"
    script:
        """
        run.py \
        --sample_id ${sample} --input_bam_path ${bam} \
        --output_path ${sample}
        """

}

workflow {
    if (!params.reference) {
        exit 1, "Reference genome not specified! Please, provide --reference"
    }

    reads = Channel
    .fromFilePairs('1-Input/*R{1,2}*.fastq')

    bwa_index = file(params.reference + ".{,amb,ann,bwt,pac,sa}")
    BWAmem(reads, bwa_index)
    xhla(BWAmem.out)
}
