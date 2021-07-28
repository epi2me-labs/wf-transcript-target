# Transcript Target Workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow
to review and consolidates transcripts of interest from direct RNA sequencing collections. 
This workflow can also be used for more general assessment of transcripts of interest.


## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-transcript-target --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* a simple text file providing a summary of sequencing reads,
* an HTML report document detailing the primary findings of the workflow,
* a BAM alignment file of all input reads with reference,
* fasta file containing consensus sequence of input sequences per reference,
* a BAM bam file of final consensus aligned with reference

## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
