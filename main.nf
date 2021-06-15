#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage(){
    log.info """
Wf-transcript-target

Usage:
    nextflow run epi2melabs/wf-transcript-target [options]

Script Options:
    --fastq             DIR     FASTQ files (required)
    --reference         FILE     Reference FASTA file (required)
    --out_dir           DIR     Path for output (default: $params.out_dir)
    --prefix            STR     The prefix attached to each of the output filenames (optional)
    --threads           INT     Number of threads per process for alignment and sorting steps (4)
    --threshold        INT     Percentage expected for consensus accuracy (85)
    --bam               BOOL    If false, bam files will not be made available in output (default: false)
    --help
    
"""
}


process fastcatQuality {
    label "wftranscripttarget"
    cpus params.threads
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
    cpus params.threads
    input:
        file "reads_*.fastq"
        file reference

    output:
        path "alignment.sam", emit: alignment
        path "alignmentStats.tsv", emit: alignmentStats
        path "readsAligned.bam", emit: alignmentBam
        path "readsAligned.bam.bai", emit: indexed
        path "*REF_*.bam", emit: splitBam
        path "unmapped-per-read.txt", emit: unmapped
        
    """
    minimap2 -t $task.cpus -ax map-ont $reference *.fastq > alignment.sam
    samtools flagstat alignment.sam -O tsv > alignmentStats.tsv
    samtools sort alignment.sam -o readsAligned.bam --threads $task.cpus
    samtools index readsAligned.bam
    bamtools split -in readsAligned.bam -mapped
    mv readsAligned.MAPPED.bam alignedReads.bam
    bamtools split -in alignedReads.bam -reference
    bedtools bamtofastq -i *UNMAPPED.bam -fq unmapped.fq
    fastcat -f unmmapped-file-summary.txt -r unmapped-per-read.txt unmapped.fq
    
"""
}

process createTuples {
    label "wftranscripttarget"
    cpus params.threads
    input:
        each file(reg)
        file reference
    output:
        tuple val("$refname"), path("*.sam"), path("*.fastq"), path("*.fasta"), emit: eachAlignment
        path "*.tsv", emit: alignmentStats
        path "*.bed.gz", emit: bedFile
    script:
        refname = "$reg".split(/\./)[1].substring(4);
        
    """
    samtools flagstat $reg -O tsv > "$refname"alignmentStats.tsv
    samtools view -h -o "$refname".sam $reg
    samtools fastq $reg > "$refname".fastq
    zcat -f $reference | seqkit grep -r -p ^"$refname" > "$refname".fasta
    samtools sort $reg -o "$refname".bam --threads $task.cpus
    samtools index "$refname".bam
    mosdepth -n --fast-mode --by 5 $refname "$refname".bam
    
    
    """

}

process consensusSeq {
    label "wftranscripttarget"
    cpus params.threads
    input:
        tuple val(refname), path(alignment), path(reads), path(reference)
    output:
        path "*Consensus.fasta", emit: seq
        path "*consensusAligned.bam", emit: alignment
        path "*consensusAligned.bam.bai", emit: indexed
        path "$reference", emit: reference
        val "$refname", emit: refname
        path "*.bed.gz", emit: bedFile
    """
    racon -t $task.cpus $reads $alignment $reference > "$refname"Consensus.fasta
    minimap2 -t $task.cpus -ax map-ont $reference "$refname"Consensus.fasta > consensusAligned.sam
    samtools sort consensusAligned.sam  -o "$refname"consensusAligned.bam --threads $task.cpus
    samtools index "$refname"consensusAligned.bam
    mosdepth -n --fast-mode --by 500 $refname "$refname"consensusAligned.bam
    """
}


process assessAssembly {
    label "wftranscripttarget"
    cpus params.threads
    input:
        file consensusSequence
        file reference
        val refname
    output:
        path "*_stats.txt", emit: stats
    """
    assess_assembly -i $consensusSequence -r $reference -p "$refname" > result.txt
   
    """
}


process report {
    label "wftranscripttarget"
    cpus 1
    input:
        file "assembly_stats/*"
        file "alignment_stats/*"
        file reference
        file alignStats
        file "consensus_seq/*"
        file qualityPerRead
        file "bedFile/*"
        file unmapped
    output:
        path "wf-transcript-target.html", emit: report
    """
    report.py wf-transcript-target.html $alignStats $qualityPerRead ${params.threshold} \
    $reference --consensus consensus_seq/* \
    --revision $workflow.revision --commit $workflow.commitId \
    --summaries assembly_stats/* --flagstats alignment_stats/* \
    --bedFiles bedFile/* --unmapped $unmapped
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
            .fromPath("${fastq}{**,.}/*.fastq*", glob: true)
            .collect()

        // Check overall quality
        quality = fastcatQuality(fastq_files)

        // Align the ref and samples and output one bam per reference
        alignments = alignReads(fastq_files, reference)

        // Create tuples with data needed for Racon(name, fastq, sam, fasta)
        seperated = createTuples(alignments.splitBam, reference)

        // Find consensus for each reference
        consensus = consensusSeq(seperated.eachAlignment)

        // Assess consensus vs reference
        assemblyStats = assessAssembly(consensus.seq,
                                       consensus.reference,
                                       consensus.refname)

        // output optional bam alignment files
        consensusAlignment = null
        consensusIndex = null
        if (params.bam == false) {
            println("")
            }
        else {
            consensusAlignment = consensus.alignment
            consensusIndex = consensus.indexed
            }
        // emit results
        alignmentBam = alignments.alignmentBam
        alignmentIndex = alignments.indexed
        consensusSeq = consensus.seq
        // create report
        report = report(assemblyStats.stats.collect(),
                        seperated.alignmentStats.collect(),
                        reference,
                        alignments.alignmentStats,
                        consensus.seq.collect(),
                        quality.perRead,
                        seperated.bedFile.collect(),
                        alignments.unmapped
                        )

        results = alignmentBam.concat(
            alignments.alignmentBam,
            alignmentIndex,
            consensusSeq,
            report,
            consensusAlignment,
            consensusIndex,
            )

    emit:
        results
}

process output {
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*", saveAs: { 
        f -> params.prefix ? "${params.prefix}-${f}" : "${f}" }
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
    reference = file(params.reference, type: "dir", checkIfExists: true)

    // Run Bioinformatics pipeline
    results = pipeline(reference, fastq)

    // output files
    output(results)

}
