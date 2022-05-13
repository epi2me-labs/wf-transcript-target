#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 
include { start_ping; end_ping } from './lib/ping'

process fastcatQuality {
    label "wftranscripttarget"
    cpus params.threads
    input: 
         tuple path(directory), val(sample_id), val(type)
    output:
        path "per-read.txt", emit: perRead
        path "*.fastq", emit: concat_reads
    """
    fastcat -s ${sample_id} -r per-read.txt -x ${directory} > all_reads.fastq
    """
}


process alignReads {
    label "wftranscripttarget"
    cpus params.threads
    input:
        file concat_reads
        file reference

    output:
        path "alignment.sam", emit: alignment
        path "alignmentStats.tsv", emit: alignmentStats
        path "readsAligned.bam", emit: alignmentBam
        path "readsAligned.bam.bai", emit: indexed
        path "*REF_*.bam", emit: splitBam
        path "unmapped-per-read.txt", emit: unmapped
        
    """
    minimap2 -t $task.cpus -ax map-ont $reference $concat_reads > alignment.sam
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
        path reg
        file "reference"
    output:
        tuple val("$refname"), path("*.sam"), path("*.fastq"), path("*.fasta"), emit: eachAlignment
        path "*.tsv", emit: alignmentStats
        path "*.bed.gz", emit: bedFile
    script:
        refname = "$reg".split(/\./)[1].substring(4);
        
    """
    echo $reg
    samtools flagstat $reg -O tsv > "$refname".alignmentStats.tsv
    samtools view -h -o "$refname".sam $reg
    samtools fastq $reg > "$refname".fastq
    seqkit grep -r -p ^"$refname" $reference > "$refname".fasta
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


process getParams {
    label "wftranscripttarget"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process getVersions {
    label "wftranscripttarget"
    cpus 1
    output:
        path "versions.txt"
    
    """
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bamtools --version | head -2 | tail -1 | sed 's/ /,/' >> versions.txt
    bedtools --version | sed 's/ /,/' >> versions.txt
    seqkit version | sed 's/ /,/' >> versions.txt
    mosdepth --version | sed 's/ /,/'  >> versions.txt
    racon --version | sed 's/^/racon,/'  >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    python -c "import pomoxis; print(pomoxis.__version__)" | sed 's/^/pomoxis,/'  >> versions.txt
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
        file versions
        path "params.json"
    output:
        path "wf-transcript-target-*.html", emit: report
    script:
        report_name = "wf-transcript-target-" + params.report_name + '.html'

    """
    report.py $report_name $alignStats $qualityPerRead ${params.threshold} \
    $reference --consensus consensus_seq/* \
    --revision $workflow.revision --commit $workflow.commitId \
    --summaries assembly_stats/* --flagstats alignment_stats/* \
    --bedFiles bedFile/* --unmapped $unmapped \
    --versions $versions \
    --parameters params.json
    """
}


// workflow module
workflow pipeline {
    take:
        reference
        fastq
    main:

        // Check overall quality
        quality = fastcatQuality(fastq)
 
        // Align the ref and samples and output one bam per reference
        alignments = alignReads(quality.concat_reads, reference)
   

        // Create tuples with data needed for Racon(name, fastq, sam, fasta)
        seperated = createTuples(alignments.splitBam.flatten(), reference)

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
        // software versions
        software_versions = getVersions()
        params_json = getParams()
        // create report
        report = report(assemblyStats.stats.collect(),
                        seperated.alignmentStats.collect(),
                        reference,
                        alignments.alignmentStats,
                        consensus.seq.collect(),
                        quality.perRead,
                        seperated.bedFile.collect(),
                        alignments.unmapped,
                        software_versions,
                        params_json
                        )

        results = alignmentBam.concat(
            alignmentIndex,
            consensusSeq,
            report,
            consensusAlignment,
            consensusIndex,
            )

    emit:
        results
        telemetry = params_json
}



process output {
    label "wftranscripttarget"
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
WorkflowMain.initialise(workflow, params, log)
workflow {
    // Start ping
    start_ping()

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

    fastq = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet, params.sanitize_fastq)

    // Acquire reference file
    reference = file(params.reference, type: "dir", checkIfExists: true)

    // Run Bioinformatics pipeline
    results = pipeline(reference, fastq)

    // output files
    output(results[0])

    // End ping
    end_ping(pipeline.out.telemetry)

}
