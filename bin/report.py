#!/usr/bin/env python
"""Create workflow report."""

import argparse

from aplanat import report
from aplanat.components import fastcat
import aplanat.graphics
from bokeh.layouts import layout
import conda_versions
import numpy as np
import pandas as pd


def read_files(summaries, sep='\t'):
    """Read a set of files and join to single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep=sep))
    return pd.concat(dfs)


def read_flag_stat_files(flagstats, sep='\t'):
    """Read a set of flag stat files and extract total alignments."""
    flagStatList = []
    for fname in sorted(flagstats):
        alignStats = pd.read_csv(fname, sep='\t', header=None)
        alignStats.reset_index()
        totalSeq = int(alignStats.iat[0, 0])
        flagStatList.append(totalSeq)
    return flagStatList


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
    tableConsensus = {'Reference name': refNames,
                      'Consensus Accuracy %': accuracyList,
                      'Number of reads aligned': flag_stats,
                      'Total Aligned %': percentageAligned}
    consensus_df = pd.DataFrame(tableConsensus)
    consensus_df[''] = np.where(
        consensus_df['Consensus Accuracy %'] < 98.0,
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
    # Assembly
    section = report_doc.add_section()
    section.markdown("## Assembly stats")
    section.markdown(
        "The following summarises the statistics from the consensus assembly"
        "with the reference")
    section.table(
        seq_summary, index=False)
    section = report_doc.add_section()

    section.markdown("## Software versions")
    section.markdown('''The table below highlights versions
                    of key software used within the analysis''')
    req = [
        'minimap2', 'samtools', 'racon', 'pomoxis', 'fastcat', 'bamtools']
    versions = conda_versions.scrape_data(
        as_dataframe=True, include=req)
    section.table(versions[['Name', 'Version', 'Build']], index=False)
    section = report_doc.add_section()

    section.markdown('''
### About

**Oxford Nanopore Technologies products are not intended for use for health
assessment or to diagnose, treat, mitigate, cure or prevent any disease or
condition.**
---
''')
    # write report
    report_doc.write(args.output)


if __name__ == "__main__":
    main()
