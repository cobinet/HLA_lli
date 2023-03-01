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

process MAKE_REF_IDS {
    container 'docker://zlskidmore/hla-la:latest'
    input:
        tuple val(sample), file(bam)
    output:
        path "ref_id_unix.txt"
    script:
        ref = "ref_id_unix.txt"
        """
        echo "contigID\tcontigLength\tExtractCompleteContig\tPartialExtraction_Start\tPartialExtraction_Stop" > ref_id.txt
        samtools idxstats ${bam} >> ref_id.txt
        awk '{ sub("\\r\$", ""); print }' ref_id.txt > ${ref}

        if grep -q "chr6\t" "${ref}"; then
            sed -i '/^chr6\t/ s/\$/28510120\t33480577/' ${ref}
        else
            sed -i '/^*\t/ s/\$/0\t0/' ${ref}
        fi
        """
}

process TYPING {
    cpus 20
    memory '200G'
    publishDir "3-HLA-LA", mode: 'link'

    graph_path = "/usr/local/bin/HLA-LA/graphs/"
    container 'docker://zlskidmore/hla-la:latest'
    containerOptions "-B ${graph}:${graph_path}/${graph} -B ${ref_ids}:${graph_path}/${graph}/knownReferences/${ref_ids}"

    input:
        tuple val(sample), file(bam)
        path graph
        path ref_ids
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
        GET_GRAPH()
        TYPING(bam, GET_GRAPH.out, MAKE_REF_IDS(GET_GRAPH.out))
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
