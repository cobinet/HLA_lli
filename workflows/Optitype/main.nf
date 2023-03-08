process OPTITYPE {
    tag "$sample"

    container 'fred2/optitype'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${prefix}"), emit: output

    script:
        prefix    = task.ext.prefix ?: "${sample}"

        """
        cat <<-END_CONFIG > config.ini
        [mapping]
        razers3=razers3
        threads=$task.cpus
        [ilp]
        threads=1
        [behavior]
        deletebam=true
        unpaired_weight=0
        use_discordant=false
        END_CONFIG

        # Run OptiType
        OptiTypePipeline.py -i ${bam} -c config.ini --${meta.seq_type} --prefix $sample --outdir $sample
        """
}