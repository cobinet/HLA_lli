include { BWA_MEM } from '../alignment/bwa_alignment.nf'

src_path = "/usr/local/bin/HLA-LA"
params.imgt_url = "https://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz"
params.hlala_ref_url = "https://dl.dropboxusercontent.com/s/mnkig0fhaym43m0/reference_HLA_ASM.tar.gz?dl=1"

process PREPARE_GRAPH {
    cpus 8
    memory '20G'

    beforeScript "mkdir ${graph.simpleName} ${reference.simpleName}"

    container 'docker://zlskidmore/hla-la:latest'
    publishDir "graph"
    containerOptions """-B ${graph.simpleName}:${src_path}/graphs/${graph.simpleName}\
                        -B ${reference.simpleName}:${src_path}/src/${reference.simpleName}"""
    input:
        path graph
        path reference
    output:
        path "$graph.simpleName", emit: graph
        path "$reference.simpleName", emit: reference
    script:
        """
        tar -xvzf $graph
        tar -xvzf $reference
        HLA-LA --action prepareGraph --PRG_graph_dir ${graph.simpleName}
        """
}

process MAKE_REF_IDS {
    container 'docker://zlskidmore/hla-la:latest'
    input:
        tuple val(sample), path(files)
    output:
        path "ref_id_unix.txt"
    script:
        ref = "ref_id_unix.txt"
        def (bam, bai) = files
        """
        echo "contigID\tcontigLength\tExtractCompleteContig\tPartialExtraction_Start\tPartialExtraction_Stop" > ref_id.txt
        samtools idxstats ${bam} | awk -v OFS='\\t' '{print \$1, \$2, 0}' >> ref_id.txt
        # Fix line returns for unix
        awk '{ sub("\\r\$", ""); print }' ref_id.txt > ref_id.unix.txt
        # Select reads to be extracted from bam
        awk -v OFS='\\t' '/^chr6_/ || /\\*/{print \$1, \$2, 1; next} 1' ref_id.unix.txt | sed '/^chr6\t/ s/\$/\t28510120\t33480577/' > ${ref}
        """
}

process TYPING {
    tag "$sample"
    cpus 20
    memory '50G'
    publishDir "3-HLA-LA"

    container 'docker://zlskidmore/hla-la:latest'
    containerOptions """-B ${graph}:${src_path}/graphs/${graph}\
                        -B ${ref_ids}:${src_path}/graphs/${graph}/knownReferences/${ref_ids}\
                        -B ${reference}:${src_path}/src/${reference}"""

    input:
        tuple val(sample), path(files)
        path graph
        path reference
        path ref_ids
    output:
        path "./out/$sample/"
    script:
        // Only alphanumeric values are allowed
        def id = (sample =~ /(\w+)/)[0][1]
        def (bam, bai) = files
        """
        mkdir out
        HLA-LA.pl\
            --BAM ${bam}\
            --graph ${graph}\
            --sampleID ${id}\
            --maxThreads ${task.cpus}\
            --workingDir ./out/\
        """
}

workflow HLA_LA {
    take: bam
    main:
        graph = Channel.fromPath(params.imgt_url)
        ref = Channel.fromPath(params.hlala_ref_url)
        PREPARE_GRAPH(graph, ref)
        TYPING(bam, PREPARE_GRAPH.out.graph, PREPARE_GRAPH.out.reference, MAKE_REF_IDS(bam))
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
