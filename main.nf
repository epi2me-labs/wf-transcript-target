#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage(){
    log.info """
Wf-transcript-target

Usage:
    nextflow run epi2melabs/wf-transcript-target [options]

Script Options:
    --fastq             DIR     FASTQ file (required)
    --reference         FILE    Reference FASTA file (required)
    --out_dir           DIR     Path for output (default: $params.out_dir)
    --bam               BOOL    If false, bam files will not be made available in output (default: false)
    --help
    
"""
}

process fastcatQuality {
    // publish inputs to output directory
    label "wftranscripttarget"
    input:
        file "reads_*.fastq"
    output:
        path "per-read.txt", emit: perRead
    """
    fastcat -f file-summary.txt -r per-read.txt *.fastq 
    """
}

process alignReads {
    label "wftranscripttarget"
    cpus 1
    input:
        file "reads_*.fastq"
        file reference

    output:
        path "alignment.sam", emit: alignment
        path "alignmentStats.tsv", emit: alignmentStats
        path "readsAligned.bam", emit: alignmentBam
        path "readsAligned.bam.bai", emit: indexed
    """
    minimap2 -ax map-ont $reference *.fastq > alignment.sam 
    samtools flagstat alignment.sam -O tsv > alignmentStats.tsv
    samtools sort alignment.sam -o readsAligned.bam
    samtools index readsAligned.bam


"""   
}

process createConsensus {
    label "wftranscripttarget"
    cpus 1
    input:
        file alignment
        file "reads_*.fastq"
        file reference
       
    output:
        path "consensus.fasta", emit: seq
        path "consensusAligned.bam", emit: alignment
        path "consensusAligned.bam.bai", emit: indexed
    """
    cat reads_*fastq > reads.fastq 
    racon reads.fastq $alignment $reference > consensus.fasta 
    minimap2 -ax map-ont $reference consensus.fasta > consensusAligned.sam 
    samtools sort consensusAligned.sam  -o consensusAligned.bam
    samtools index consensusAligned.bam
    """
}

process assessAssembly {
    label "wftranscripttarget"
    cpus 1
    input:
        file consensusSequence
        file reference
       
    output:
        path "assemblyResult_stats.txt", emit: stats
    """
    assess_assembly -i $consensusSequence -r $reference -p assemblyResult > result.txt 

    """
}


process report {
    label "wftranscripttarget"
    cpus 1
    input:
        file assemblyStats
        file alignStats
        file qualityPerRead
    output:
        path "wf-transcript-target.html", emit: report
    """
    report.py wf-transcript-target.html $assemblyStats $alignStats $qualityPerRead

    """
}
// workflow module
workflow pipeline {
    take:
        reference
        fastq
    main:
        // Get fastq files from dir path
        fastq_files = channel
            .fromPath("${fastq}{**,.}/*.fastq", glob: true)
            .collect()
         //quality step
        quality = fastcatQuality(fastq_files)
        // Align the two input files and create stats
        alignments = alignReads(fastq_files, reference)

        // Using output of alignReads find consensus
        consensus = createConsensus(alignments.alignment, fastq_files, reference)
        
        // Assess consensus vs reference
        assemblyStats = assessAssembly(consensus.seq, reference)

        // report
        report = report(assemblyStats.stats,
                        alignments.alignmentStats,
                        quality.perRead)
        // output optional bam alignment files
        consensusAlignment = null
        consensusIndex = null
        if (params.bam == false) {
            println("")
        } else {
            consensusAlignment = consensus.alignment
            consensusIndex = consensus.indexed   
        }

        //emit results 
        alignmentBam = alignments.alignmentBam
        alignmentIndex = alignments.indexed
        consensusSeq = consensus.seq
    
        results = alignmentBam.concat(
            alignmentIndex,
            consensusSeq, 
            consensusAlignment,
            consensusIndex,
            report)
           

    emit:
        results
        
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wftranscripttarget"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}

// entrypoint workflow
workflow {

    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq) {
        helpMessage()
        println("")
        println("`--fastq` is required")
        exit 1
    }

    if (!params.reference) {
        helpMessage()
        println("")
        println("`--reference` is required")
        exit 1
    }

    // Acquire fastq test directory
    fastq = file(params.fastq, type: "dir", checkIfExists: true)

    // Acquire reference file
    reference = file(params.reference, type: "file", checkIfExists: true)

    // Run Bioinformatics pipeline
    results = pipeline(reference, fastq)

    // output files
    output(results)
  
}
