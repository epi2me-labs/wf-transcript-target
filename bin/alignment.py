#!/usr/bin/env python
"""Alignment script."""

from itertools import count, takewhile

from Bio import SeqIO
import edlib


def sliced(seq, n):
    """Slice up a string."""
    iterator = takewhile(len, (seq[i: i + n] for i in count(0, n)))
    return iterator


def alignment(query, target, splitnumber):
    """Alignment function Needle is default."""
    aligned = edlib.align(query, target, task="path")
    query_len = len(query)
    query_target = len(target)
    nice_aligned = edlib.getNiceAlignment(aligned, query, target)
    query = list(sliced(nice_aligned['query_aligned'], splitnumber))
    matched = list(sliced(nice_aligned['matched_aligned'], splitnumber))
    target = list(sliced(nice_aligned['target_aligned'], splitnumber))
    formatted = str('Edit Distance: '+str(aligned['editDistance'])+'\n')
    total = ''
    for i in range(0, len(query)):
        total += 'Consensus  ' + query[i] + '<br>'
        total += 'Alignment  ' + matched[i] + '<br>'
        total += 'Reference  ' + target[i] + '<br>' + '<br>'
    return total, formatted, query_len, query_target


def referenceSeq(filenames):
    """Get sequence dictionary from reference fasta file."""
    theFile = open(str(filenames))
    refDic = {}
    for seq_record in (SeqIO.parse(theFile, "fasta")):
        refDic[str(seq_record.name)] = str((seq_record.seq)).replace("T", "U")
    return refDic
