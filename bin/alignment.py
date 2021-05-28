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
    queryLen = len(query)
    queryTar = len(target)
    niceAligned = edlib.getNiceAlignment(aligned, query, target)
    query = list(sliced(niceAligned['query_aligned'], splitnumber))
    matched = list(sliced(niceAligned['matched_aligned'], splitnumber))
    target = list(sliced(niceAligned['target_aligned'], splitnumber))
    formatted = str('Edit Distance: '+str(aligned['editDistance'])+'\n')
    total = ''
    for i in range(0, len(query)):
        total += 'Consensus  ' + query[i] + '<br>'
        total += 'Alignment  ' + matched[i] + '<br>'
        total += 'Reference  ' + target[i] + '<br>' + '<br>'
    return total, formatted, queryLen, queryTar


def referenceSeq(filenames):
    """Get sequence dictionary from reference fasta file."""
    theFile = open(str(filenames))
    refDic = {}
    for seq_record in (SeqIO.parse(theFile, "fasta")):
        refDic[str(seq_record.name)] = str((seq_record.seq)).replace("T", "U")
    return refDic
