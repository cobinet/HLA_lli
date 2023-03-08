process PHLAT {
    container "docker://chrisamiller/phlat:1.1_withindex"
    input:
        tuple val(sample), path(reads)
    output:
        path "out/"

    script:
        def f1, f2 = reads
        """
        mkdir out
        python -O /opt/dist/PHLAT.py\
            -1 $f1\
            -2 $f2\
            -index /opt/b2folder\
            -b2url `which bowtie2`\
            -orientation "--fr"\
            -tag $sample\
            -e /opt\
            -o out\
        """
}