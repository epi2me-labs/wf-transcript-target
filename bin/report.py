#!/usr/bin/env python
"""Create workflow report."""

import argparse

import aplanat.graphics
from aplanat.report import HTMLReport
import pandas as pd


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument("summaries", nargs='+', help="Read summary file.")
    parser.add_argument("alignmentStats", nargs='+',
                        help="Alignment stats file.")
    args = parser.parse_args()
    exec_summary = aplanat.graphics.InfoGraphItems()
    stats = args.summaries
    statsdf = pd.read_csv(stats[0], sep='\t')
    statsdf = statsdf.drop(['rstart', 'rend'], 1)
    consensusAccuracy = statsdf.iat[0, -1]
    report = HTMLReport(
        "Workflow Transcript Target report",
        ("Results generated through the wf-transcript target nextflow"
         "workflow by Oxford Nanopore Technologies"))
    report.markdown("## Assembly stats")
    report.markdown(
        "The following summarises the statistics from the consensus assembly"
        "with the reference")
    report.table(
        statsdf, index=False)
    alignStats = args.alignmentStats
    alignStats = pd.read_csv(alignStats[0], sep='\t', header=None)
    alignStats.reset_index()
    totalSeq = int(alignStats.iat[0, 0])
    mapped = int(alignStats.iat[4, 0])
    percentMapped = (mapped/totalSeq)*100
    report.markdown("## Executive summary", key="exec-head")
    report.markdown(
        "The following summarises the key findings of this workflow.",
        key="exec-desc")
    exec_summary.append("Consensus Accuracy",
                        str("%.2f" % round(consensusAccuracy, 2)) + '%',
                        "bullseye")
    exec_summary.append("On target reads",
                        str("%.2f" % round(percentMapped, 2)) + '%', "percent")
    report.plot(None, "exec-plot")
    exec_plot = aplanat.graphics.infographic(exec_summary.values(), ncols=4)
    report.plot(exec_plot, key="exec-plot")

    report.markdown('''
### About

**Oxford Nanopore Technologies products are not intended for use for health
assessment or to diagnose, treat, mitigate, cure or prevent any disease or
condition.**

This report was produced using the
[epi2me-labs/wf-template](https://github.com/epi2me-labs/wf-template).  The
workflow can be run using `nextflow epi2me-labs/wf-template --help`
---
''')
    # write report
    report.write(args.report)


if __name__ == "__main__":
    main()
