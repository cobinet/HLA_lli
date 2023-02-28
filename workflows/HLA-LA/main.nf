include { BWA_MEM } from '../alignment/bwa_alignment.nf'

process GET_GRAPH {
    cpus 20
    memory '20G'
    container 'docker://zlskidmore/hla-la:latest'
    publishDir "graph"
    output:
        path "PRG_MHC_GRCh38_withIMGT/", emit: graph
    script:
        """
        wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz
        tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz
        HLA-LA --action prepareGraph --PRG_graph_dir PRG_MHC_GRCh38_withIMGT
        """
}

process TYPING {
    cpus 20
    memory '200G'
    container 'docker://zlskidmore/hla-la:latest'
    publishDir "3-HLA-LA", mode: 'link'
    input:
        tuple val(sample), file(bam)
        path graph
    output:
        path "./out/*"
    script:
        """
        mkdir out
        HLA-LA.pl\
            --BAM ${bam}\
            --graph ${graph}\
            --sampleID ${sample}\
            --maxThreads ${task.cpus}\
            --workingDir ./out/
        """
}

workflow HLA_LA {
    take: bam
    main:
        TYPING(bam, GET_GRAPH())
    emit:
        TYPING.out
}

workflow {
    if (!params.ref) {
        exit 1, "Reference genome not specified! Please, provide --ref"
    }

    reads = Channel
    .fromFilePairs('1-Input/*R{1,2}*.fastq')

    bwa_index = file(params.ref + "*")
    BWA_MEM(reads, bwa_index)
    HLA_LA(BWAmem.out)
}
