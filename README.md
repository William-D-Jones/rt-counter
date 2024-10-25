# Transcriptional Read-through Counter (rt-counter)

The Transcription Read-through Counter (rt-counter) counts RNA-seq reads from
paired-end data that show evidence of transcriptional read-through, that is,
when RNA polymerase reads past the typical transcription termination site
and into the flanking genomic DNA.

This program reads features from a .gtf file containing the information about
each gene or transcript, including the site of transcriptional termination.
The program also takes as input a .bam file reporting paired-end alignments
from an RNA-seq project. Reads are counted for each gene and are classified as
either 'base' reads--any read that maps a gene--or 'readthrough' reads--base
reads that also extend past the 3' end of a gene and into the surrounding
genomic DNA.

## Detailed Read Counting Procedure:

    Counts paired-end reads from a BAM file using the following procedure:
    1. A template counts toward the 'base' transcript of a feature F if:
        a. It is mapped uniquely.
        b. It fulfills these "check" criteria according to
           the aligner: all segments are aligned, properly paired, not
           PCR or optical duplicates, and do not fail the platform QC.
        c. At least one segment in the template has at least fAnchor CIGAR
           match operations (M, =, or X) that overlap F.
        d. No CIGAR match operations in either segment of the template
           overlap any features other than F.
        e. With respect to the direction of transcription of F, no 
           CIGAR match operations in either segment of the template
           overlap any interval lying upstream of the most upstream exon of F.
        f. With respect to the direction of transcription of F, no 
           CIGAR match operations in either segment of the template
           overlap any interval lying downstream of the most upstream exon of 
           F, unless the template also fulfills condition g.
    2. A template also counts toward the 'readthrough' transcript of F if
       it counts toward F's 'base' transcript and also:
        g. At least one segment in the template has at least 1 CIGAR
           match operation that overlaps an interval lying downstream of the 
           most downstream exon of F, as long as there no skips between
           this operation and the more downstream of either: the most 
           upstream match operation in the segment or a non-skip operation
           containing the most downstream nucleotide of the most downstream
           exon in F.

        For every template satisfying criteria a-g, we calculate the
        'tail length' for each segment:
           The tail length is the length of the interval bounded
           by, but not containing, the most downstream nucleotide of the most 
           downstream exon in F and the most upstream N operation in the
           segment, with respect to F, that is downstream of the most
           downstream match operation in the segment that overlaps F. In
           other words, the interval used for the tail length must contain
           no N operations.

    In addition, if pathToHits is specified, all matching reads are written
    to this path and file. Alignments are written in SAM format with these 
    added optional fields:
    1. fe:Z:FEATURE_NAME
    2. ba:i:LEN (count of match operations supporting base transcript)
    3. tl:i:TAIL_LENGTH (this field is ommitted if there is no readthrough)

    Finally, if pathToFails is specified, all nonmatching reads are written
    to this path and file. Alignments are written in SAM format with 1 added 
    optional field:
    1. xx:i:ERROR_CODE

## Prerequisites:

Note: Other versions of Python packages may work as well.

1. Python 3.9.18

2. argparse 1.4.0

3. numpy 2.1.2

4. pandas 2.2.3

5. HTSeq 2.0.0

## Getting Started:

```
git clone https://github.com/William-D-Jones/rt-counter.git
```

## How to Count Readthrough Transcripts:

1. Before beginning, paired-end RNA-seq reads must be sorted by name so that
mates appear in adjacent records in the .bam file (in other words, reads
should not be sorted by coordinate). For example, `samtools sort -n` produces
an acceptable .bam file.

2. Use an input .gtf containing all the features for which transcriptional
read-through should be evaluated. Only features of the type "exon" are
considered, and exons are grouped by the "gene_id" field.

3. Run `rt-counter` as described below.

## Program Usage:

```
usage: rt_counter.py [-h] [-m PATH_TO_HITS] [-x PATH_TO_FAILS]
                     [-f NUMBER_OF_BASES]
                     PATH_TO_GTF PATH_TO_BAM PATH_TO_OUTPUT_COUNTS
                     PATH_TO_OUTPUT_TAILS

Counts templates showing evidence of transcriptional read-through in paired-
end RNA-seq data.

positional arguments:
  PATH_TO_GTF           Path to Ensembl GTF file containing the genomic
                        features, or any feature file readable by
                        HTSeq.GFF_Reader with end_included=True.
  PATH_TO_BAM           Path to BAM file containing the aligned reads. BAM
                        alignments must be sorted by read name (QNAME).
  PATH_TO_OUTPUT_COUNTS
                        Path and filename to which the output tsv file of
                        counting results will be written.
  PATH_TO_OUTPUT_TAILS  Path and filename to which the output tsv file of
                        readthrough tail lengths will be written.

options:
  -h, --help            show this help message and exit
  -m PATH_TO_HITS, --pathToHits PATH_TO_HITS
                        If present, writes alignments from the BAM file (no
                        headers) to this path and filename if a match is found
                        to a GTF entry, with 3 added optional fields:
                        fe:Z:FEATURE_NAME, ba:i:tail_length_base,
                        rt:i:tail_length_readthrough
  -x PATH_TO_FAILS, --pathToFails PATH_TO_FAILS
                        If present, writes alignments from the BAM file (no
                        headers) to this path and filename if the alignment is
                        not counted toward any GTF entry, with 1 added option
                        field: xx:i:ERROR_CODE
  -f NUMBER_OF_BASES, --featureAnchor NUMBER_OF_BASES
                        Number of bases in reads that must match the 3' end of
                        each feature in the GTF file. Defaults to 10.

```
