#!/usr/bin/env nextflow

// Developer notes
// 
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion.
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

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
    --help

    
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
    """
    minimap2 -ax map-ont $reference *.fastq > alignment.sam 

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
        path "consensus.fasta", emit: consensus
    """
    cat reads_*fastq > reads.fastq 
    racon reads.fastq $alignment $reference > consensus.fasta 
    """
}

process assessAssembly {
    label "wftranscripttarget"
    cpus 1
    input:
        file consensus
        file reference
       
    output:
        path "assemblyResult_stats.txt", emit: stats
        path "assemblyResult_summ.txt", emit: assemblySummary
    """
    assess_assembly -r $consensus -i $reference -p assemblyResult > result.txt 

    """
}


process report {
    label "wftranscripttarget"
    cpus 1
    input:
        file assemblyStats

    output:
        path "wf-transcript-target.html", emit: report
    """
    report.py wf-transcript-target.html $assemblyStats 

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
        // Align the two input files and create stats
        alignment = alignReads(fastq_files, reference)

        //Using output of alignReads find consensus
        consensus = createConsensus(alignment,fastq_files,reference)

        //Assess consensus vs reference
        assemblyStats = assessAssembly(consensus,reference)

        //report
        report = report(assemblyStats.stats)
        
    emit:
        alignment = alignment
        consensus = consensus
        stats = assemblyStats.stats
        summmary = assemblyStats.assemblySummary
        report = report
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
    reference = file(params.reference, type:"file", checkIfExists:true)

    // Run Bioinformatics pipeline
    results = pipeline(reference,fastq)

    //output(results.consensus,results.alignment,results.assemblyStats.stats,results.a.assemblySummary)
    output(results.alignment.concat(
        results.consensus,results.stats,results.summmary,results.report))
}
