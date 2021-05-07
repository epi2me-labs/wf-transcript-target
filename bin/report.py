"""Create workflow report."""

import argparse

from aplanat.components import fastcat
import aplanat.graphics

from aplanat.report import HTMLReport

from bokeh.layouts import layout

import conda_versions
import pandas as pd


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file.")
    parser.add_argument("summaries", nargs='+', help="Read summary file.")
    parser.add_argument("alignmentStats", nargs='+',
                        help="Alignment stats file.")
    parser.add_argument("quality", nargs='+', help="Fastcat quality results.")
    report = HTMLReport(
        "Workflow Transcript Target report",
        ("Results generated through the wf-transcript target nextflow"
         " workflow by Oxford Nanopore Technologies"))
    args = parser.parse_args()
    # Quality
    quality = args.quality
    quality_df = pd.read_csv(quality[0], sep='\t')
    read_qual = fastcat.read_quality_plot(quality_df)
    read_length = fastcat.read_length_plot(quality_df)
    report.markdown("## Read Quality Control")
    report.markdown("This sections displays basic QC"
                    "metrics indicating read data quality.")
    report.plot(
        layout(
            [[read_length, read_qual]],
            sizing_mode="stretch_width"))

    # Assembly
    stats = args.summaries
    statsdf = pd.read_csv(stats[0], sep='\t')
    statsdf = statsdf.drop(['rstart', 'rend'], 1)
    consensusAccuracy = statsdf.iat[0, -1]
    report.markdown("## Assembly stats")
    report.markdown(
        "The following summarises the statistics from the consensus assembly"
        "with the reference")
    report.table(
        statsdf, index=False)

    # Align stats
    alignStats = args.alignmentStats
    alignStats = pd.read_csv(alignStats[0], sep='\t', header=None)
    alignStats.reset_index()
    totalSeq = int(alignStats.iat[0, 0])
    mapped = int(alignStats.iat[4, 0])
    percentMapped = (mapped/totalSeq)*100
    # Exec Summary infographic
    report.markdown("## Executive summary", key="exec-head")
    report.markdown(
        "The following summarises the key findings of this workflow.",
        key="exec-desc")
    exec_summary = aplanat.graphics.InfoGraphItems()
    exec_summary.append("Consensus Accuracy",
                        str("%.2f" % round(consensusAccuracy, 2)) + '%',
                        "bullseye")
    exec_summary.append("On target reads",
                        str("%.2f" % round(percentMapped, 2)) + '%', "percent")
    report.plot(None, "exec-plot")
    exec_plot = aplanat.graphics.infographic(exec_summary.values(), ncols=4)
    report.plot(exec_plot, key="exec-plot")
    report.markdown("## Software versions", key="version-head")
    report.markdown('''The table below highlights versions
                    of key software used within the analysis''',
                    key="version-desc")
    req = [
        'minimap2', 'samtools', 'racon', 'pomoxis', 'fastcat']
    versions = conda_versions.scrape_data(
        as_dataframe=True, include=req)
    report.table(versions[['Name', 'Version', 'Build']], index=False)
    report.markdown('''
### About

**Oxford Nanopore Technologies products are not intended for use for health
assessment or to diagnose, treat, mitigate, cure or prevent any disease or
condition.**
---
''')
    # write report
    report.write(args.report)


if __name__ == "__main__":
    main()
