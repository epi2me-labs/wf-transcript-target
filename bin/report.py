#!/usr/bin/env python
"""Create workflow report."""

import argparse

import alignment
from aplanat import report
from aplanat.components import fastcat
import aplanat.graphics
from aplanat.lines import steps
from bokeh.layouts import gridplot, layout
import conda_versions
import numpy as np
import pandas as pd


def read_files(summaries, sep='\t'):
    """Read a set of files and join to single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep=sep))
    return pd.concat(dfs)


def read_consensus(consensus, sep='\t'):
    """Read a set of conensus files and extract sequences."""
    consensus_dic = {}
    for fname in sorted(consensus):
        new_dic = alignment.referenceSeq(str(fname))
        consensus_dic.update(new_dic)
    return consensus_dic


def read_flag_stat_files(flagstats, sep='\t'):
    """Read a set of flag stat files and extract total alignments."""
    flagStatList = []
    for fname in sorted(flagstats):
        alignStats = pd.read_csv(fname, sep='\t', header=None)
        alignStats.reset_index()
        totalSeq = int(alignStats.iat[0, 0])
        flagStatList.append(totalSeq)
    return flagStatList


def depth_graph(bedFiles, sep='\t'):
    """Create depth vs position graph."""
    graphs = []
    for fname in sorted(bedFiles):
        binned_depth = pd.read_csv(
            str(fname), sep=sep,
            names=['ref', 'start', 'end', 'depth'])
        p = steps(
            [binned_depth['start']],
            [binned_depth['depth']],
            colors=['darkolivegreen'],
            x_axis_label='Position along reference',
            y_axis_label='Sequencing depth / Bases',
            title=str(fname[8:-15:]))
        p.title.align = "left"
        p.title.text_font_size = "12px"
        p.xaxis.formatter.use_scientific = False
        graphs.append(p)
    return graphs


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "output",
        help="Report output file.")
    parser.add_argument(
        "alignmentStats",
        help="Alignment stats file.")
    parser.add_argument(
        "quality",
        help="Fastcat quality results.")
    parser.add_argument(
        "threshold",
        help="Threashold percentage expected for consensus accuracy")
    parser.add_argument(
        "references",
        help="Reference file input at start")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--summaries", nargs='+',
        help="Read summary file.")
    parser.add_argument(
        "--flagstats", nargs='+',
        help="Flag stat summaries")
    parser.add_argument(
        "--consensus", nargs='+',
        help="consensus fasta sequences")
    parser.add_argument(
        "--bedFiles", nargs='+',
        help="bed files for sequence depth")
    args = parser.parse_args()

    report_doc = report.WFReport(
        "Transcript target report", "wf-transcript-target",
        revision=args.revision, commit=args.commit)

    flag_stats = read_flag_stat_files(args.flagstats)

    # Retrieve flag and consensus stats and create df
    seq_summary = read_files(args.summaries)
    statsdf = seq_summary
    statsdf = statsdf.drop(['rstart', 'rend'], 1)
    refNames, accuracyList = list(statsdf['ref']), list(statsdf['acc'])
    accuracyList = list(np.around(np.array(accuracyList), 2))
    alignStats = args.alignmentStats
    alignStats = pd.read_csv(alignStats, sep='\t', header=None)
    alignStats.reset_index()
    # Retrieve total flag stats
    totalSeq = int(alignStats.iat[0, 0])
    mapped = int(alignStats.iat[4, 0])
    percentMapped = (mapped/totalSeq)*100
    percentageAligned = list(map((lambda x: (x/totalSeq) * 100), flag_stats))
    percentageAligned = list(np.around(np.array(percentageAligned), 2))
    # Output all in a table
    threshold = int(args.threshold)
    tableConsensus = {'Reference name': refNames,
                      'Consensus Accuracy %': accuracyList,
                      'Number of reads aligned': flag_stats,
                      'Total Aligned %': percentageAligned}
    consensus_df = pd.DataFrame(tableConsensus)
    consensus_df[''] = np.where(
        consensus_df['Consensus Accuracy %'] < threshold,
        'Warning', '')
    section = report_doc.add_section()
    section.markdown("## Summary")
    section.markdown(" This table summarises the consensus accuracy",
                     "and alignment stats for each reference.")
    section.table(consensus_df, index=False)
    if consensus_df[''].isin({'': ['Warning']}).any():
        section.markdown('**Warning: Some references < threshold accuracy**')
    else:
        pass
    # Exec Summary infographic
    section = report_doc.add_section()
    section.markdown("## Total aligned reads")
    otr = [("On target reads",
            str("%.2f" % round(percentMapped, 2)) + '%',
            "percent", '')]
    exec_plot = aplanat.graphics.infographic(otr, ncols=1)
    section.plot(exec_plot, key="exec-plot")
    # Quality
    quality = args.quality
    quality_df = pd.read_csv(quality, sep='\t')
    read_qual = fastcat.read_quality_plot(quality_df)
    read_length = fastcat.read_length_plot(quality_df)
    section = report_doc.add_section()
    section.markdown("## Read Quality Control")
    section.markdown("This sections displays basic QC"
                     "metrics indicating read data quality.")
    section.plot(
        layout(
            [[read_length, read_qual]],
            sizing_mode="stretch_width"))
    # Depth Coverage
    section = report_doc.add_section()
    section.markdown("## Depth of coverage")
    section.markdown("The depth of coverage of alignments "
                     "across the reference.")
    depth = args.bedFiles
    depthGraphs = depth_graph(depth)
    section.plot(
            gridplot(depthGraphs, ncols=2, plot_width=300, plot_height=300))
    # Assembly
    section = report_doc.add_section()
    section.markdown("## Consensus Alignment Statistics")
    section.markdown(
        "The following summarises the statistics from the consensus"
        "aligned with the reference")
    section.table(
        seq_summary, index=False)
    # Reference and consensus alignments
    section = report_doc.add_section()
    section.markdown("## Alignment of Consensus and Reference")
    section.markdown(
        "Sequence alignment using Levenshtein (edit) distance.")
    refFile = args.references
    refseq = alignment.referenceSeq(str(refFile))
    cons = read_consensus(args.consensus)
    for item in (cons.keys()):
        section.markdown("###" + str(item))
        p = alignment.alignment(cons[str(item)], refseq[str(item)], 50)
        section.markdown(str(p[1]))
        section.markdown('Length of Consensus: ' + str(p[2]))
        section.markdown('Length of Reference: ' + str(p[3]))
        section.markdown("<pre>" + p[0] + "</pre>")
    section = report_doc.add_section()
    section.markdown("## Software versions")
    section.markdown('''The table below highlights versions
                    of key software used within the analysis''')
    req = [
        'minimap2', 'samtools', 'racon', 'pomoxis', 'fastcat', 'bamtools',
        'python-edlib', 'biopython', 'mosdepth']
    versions = conda_versions.scrape_data(
        as_dataframe=True, include=req)
    section.table(versions[['Name', 'Version', 'Build']], index=False)
    section = report_doc.add_section()
    # write report
    report_doc.write(args.output)


if __name__ == "__main__":
    main()
