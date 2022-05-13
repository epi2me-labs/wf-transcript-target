## Introduction

This workflow is suitable for generating summary information from a small number of user-supplied transcripts
of interest. 

In order to do this, an alignment of direct RNA reads to the reference transcripts is done with [minimap2](https://github.com/lh3/minimap2).
Reads mapping to the reference transcripts are used to create consensus assemblies using
[racon](https://github.com/isovic/racon), which are then compared to the original reference sequences. 




