#!/usr/bin/env python
"""Create workflow report."""

import argparse

import alignment
from aplanat import report
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
import aplanat.graphics
from aplanat.lines import steps
from bokeh.layouts import gridplot, layout
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
    flag_stat_list = []
    for fname in sorted(flagstats):
        align_stats = pd.read_csv(fname, sep='\t', header=None)
        align_stats.reset_index()
        total_seq = int(align_stats.iat[0, 0])
        flag_stat_list.append(total_seq)
    return flag_stat_list


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
        "alignment_stats",
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
    parser.add_argument(
        "--unmapped", nargs='+',
        help="unmapped fastcat stats")
    parser.add_argument(
        "--versions", required=True,
        help="CSV containing name, version.")
    parser.add_argument(
        "--parameters", required=True,
        help="JSON of parameters used for wf")
    args = parser.parse_args()
    report_doc = report.WFReport(
        "Transcript target report", "wf-transcript-target",
        revision=args.revision, commit=args.commit)

    flag_stats = read_flag_stat_files(args.flagstats)
    # Retrieve flag and consensus stats and create df
    seq_summary = read_files(args.summaries)
    statsdf = seq_summary
    statsdf = statsdf.drop(['rstart', 'rend'], 1)
    ref_names, accuracy_list = list(statsdf['ref']), list(statsdf['acc'])
    accuracy_list = list(np.around(np.array(accuracy_list), 2))
    align_stats = args.alignment_stats
    align_stats = pd.read_csv(align_stats, sep='\t', header=None)
    align_stats.reset_index()
    # Retrieve total flag stats
    total_seq = int(align_stats.iat[0, 0])
    mapped = int(align_stats.iat[4, 0])
    percent_mapped = (mapped/total_seq)*100
    percentage_aligned = list(map((lambda x: (x/total_seq) * 100), flag_stats))
    percentage_aligned = list(np.around(np.array(percentage_aligned), 2))
    # Output all in a table
    threshold = int(args.threshold)
    table_consensus = {
        'Reference name': ref_names,
        'Consensus Accuracy %': accuracy_list,
        'Number of reads aligned': flag_stats,
        'Total Aligned %': percentage_aligned}
    consensus_df = pd.DataFrame(table_consensus)
    consensus_df[''] = np.where(
        consensus_df['Consensus Accuracy %'] < threshold,
        'Warning', '')
    section = report_doc.add_section()
    section.markdown("## Summary")
    section.markdown(
        " This table summarises the consensus accuracy",
        "and alignment stats for each reference.")
    section.table(consensus_df, index=False)
    if consensus_df[''].isin({'': ['Warning']}).any():
        section.markdown('**Warning: Some references < threshold accuracy**')
    else:
        pass
    # Exec Summary infographic
    unmapped = args.unmapped
    print(unmapped)
    unmapped_df = pd.read_csv(unmapped[0], delimiter='\t')
    unmapped_read_qual = unmapped_df["mean_quality"].mean()
    unmapped_read_length = unmapped_df["read_length"].mean()
    section = report_doc.add_section()
    section.markdown("## Total aligned reads")
    exec_summary = aplanat.graphics.InfoGraphItems()
    exec_summary.append(
        "Total reads",
        str(total_seq),
        "calculator", '')
    exec_summary.append(
        "On target reads",
        str("%.2f" % round(percent_mapped, 2)) + '%',
        "percent", '')
    exec_summary.append(
        "Unmapped mean Q",
        str("%.2f" % round(unmapped_read_qual, 2)),
        "clipboard-check", '')
    exec_summary.append(
        "Unmapped mean len",
        str("%.0f" % round(unmapped_read_length, 0)),
        "calculator", '')
    exec_plot = aplanat.graphics.infographic(exec_summary.values(), ncols=4)
    section.plot(exec_plot, key="exec-plot")
    # Quality
    quality = args.quality
    quality_df = pd.read_csv(quality, sep='\t')
    read_qual = fastcat.read_quality_plot(quality_df)
    read_length = fastcat.read_length_plot(quality_df)
    section = report_doc.add_section()
    section.markdown("## Read Quality Control")
    section.markdown(
        "This sections displays basic QC"
        " metrics indicating read data quality.")
    section.plot(
        layout(
            [[read_length, read_qual]],
            sizing_mode="stretch_width"))
    # Depth Coverage
    section = report_doc.add_section()
    section.markdown("## Depth of coverage")
    section.markdown(
        "The depth of coverage of alignments "
        "across the reference.")
    depth = args.bedFiles
    depth_graphs = depth_graph(depth)
    section.plot(
            gridplot(depth_graphs, ncols=2, plot_width=300, plot_height=300))
    # Assembly
    section = report_doc.add_section()
    section.markdown("## Consensus Alignment Statistics")
    section.markdown(
        "The following summarises the statistics from the consensus"
        "aligned with the reference")
    seq_summary = statsdf.drop(['name'], 1)
    seq_summary = seq_summary.rename(columns={'ref': 'Name'})
    seq_summary = seq_summary.round({'ref_coverage': 4, 'iden': 4, 'acc': 4})
    section.table(seq_summary)
    # Reference and consensus alignments
    section = report_doc.add_section()
    section.markdown("## Alignment of Consensus and Reference")
    section.markdown(
        "Sequence alignment using Levenshtein (edit) distance.")
    ref_file = args.references
    refseq = alignment.referenceSeq(str(ref_file))
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
    section = report_doc.add_section(
        section=scomponents.version_table(args.versions))
    # Params reporting
    report_doc.add_section(
        section=scomponents.params_table(args.parameters))
    # write report
    report_doc.write(args.output)


if __name__ == "__main__":
    main()
