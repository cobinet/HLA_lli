process PHLAT {
    container "docker://chrisamiller/phlat:1.1_withindex"
    publishDir "2-PHLAT"
    input:
        tuple val(sample), path(reads)
    output:
        path "out/*.sum"

    script:
        def (f1, f2) = reads
        """
        mkdir out
        python2 -O /opt/dist/PHLAT.py -1 ${f1}\
            -2 ${f2}\
            -index /opt/b2folder\
            -b2url `which bowtie2`\
            -orientation "--fr"\
            -tag ${sample}\
            -e /opt\
            -o out\
        """
}

workflow {
    files = Channel.fromFilePairs("test-chr6/*{1,2}.fastq")
    files.view()
    PHLAT(files) | view
}