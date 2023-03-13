process POLYSOLVER {
    container "vacation/hla-polysolver:latest"
    input:
        val sample
        path bam
        path bai
    output:
        "output/"
    script:
        race = "Caucasian"
        includeFreq = 1
        insertCalc = 0
        fromat = "STDFQ"
        ref = "hg19" // TODO: can this be hg38
        """
        mkdir output
        shell_call_hla_type\
            $bam\
            $race\
            $includeFreq\
            $ref\
            $format\
            $insertCalc\
            output
        """
}

workflow {
    bam = Channel.fromPath("alignment/*.bam")
    bai = Channel.fromPath("alignment/*.bai")
    sample = bam.map({it ->
        it.baseName()
    })
    POLYSOLVER(tuple(sample, bam, bai))
}