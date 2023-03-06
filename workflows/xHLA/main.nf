include { BWA_MEM } from '../alignment/bwa_alignment.nf'

process XHLA {
    cpus 20
    memory '500G'
    container 'docker://humanlongevity/hla'
    publishDir "3-xHLA"
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
    BWA_MEM(reads, bwa_index)
    XHLA(BWA_MEM.out)
}
