process TYPING {
    tag "$sample"
    cpus 10
    memory '20G'

    container 'quay.io/biocontainers/optitype:1.3.5--0'
    publishDir "2-Optitype", mode: 'symlink'

    input:
        tuple val(sample), path(files)

    output:
        tuple val(sample), path("out/")

    script:

        """
        cat <<-END_CONFIG > config.ini
        [mapping]
        razers3=razers3
        threads=$task.cpus
        [ilp]
        solver=glpk
        threads=1
        [behavior]
        check_sq=false
        deletebam=false
        unpaired_weight=0
        use_discordant=false
        END_CONFIG

        # Run OptiType
        OptiTypePipeline.py -i $files -c config.ini -d --prefix $sample --outdir out
        """
}

workflow OPTITYPE {
    take: reads
    main:
        TYPING(reads)
    emit:
        TYPING.out
}

workflow {
    files = Channel.fromFilePairs('1-Input/*_R{1,2}*')
    TYPING(files) | view
}