#!/usr/bin/env python
"""Create workflow report."""

import argparse

from aplanat.report import HTMLReport
import pandas as pd


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument("summaries", nargs='+', help="Read summary file.")
    args = parser.parse_args()
    stats = args.summaries
    statsdf = pd.read_csv(stats[0], sep='\t')
    statsdf = statsdf.drop(['ref', 'rstart', 'rend', 'ref_coverage'], 1)
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