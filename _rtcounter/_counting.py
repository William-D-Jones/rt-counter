import HTSeq as ht

# Parser choices
MODE_CNTALL="count_all"
MODE_CNTNON="count_none"
MODE_CNTFRA="count_fractional"
COUNTMODES=[MODE_CNTALL,MODE_CNTNON,MODE_CNTFRA]

def _calc_rt(pathToGTF,pathToBAM,pathToHits,pathToFails,fAnchor,rtAnchor,\
        countMultimappers,countOverlappers):
    """
    Counts reads from a BAM file using the following procedure:
    1. If a template contains at least fAnchor nucleotides anywhere within a
    feature, we count the template toward the feature's 'base' transcript.
    2. For each match in (1), if the template matches at least fAnchor
    nucleotides counting back from the 3' end of the feature, we count the
    template toward the feature's 'end' transcript.
    3. For each match in (2), if the template also matches at least the first
    fAnchor nucleotides counting forward from the 3' end of the feature (that
    is, the genomic DNA downstream of the feature), we count the template
    toward the feature's 'readthrough' transcript.

    The counting procedure is performed at the exon level. Therefore, if 
    multiple transcripts are produced from the same genomic feature, only
    readthrough in transcripts that contain the full length of the most 3'
    exon (relative to the direction of transcription) will have their 
    readthrough correctly evaluated.

    Returns a Pandas dataframe containing these columns, left-to-right:
    0. name: the name of the gene
    1. counts_base: the counts to the base transcript
    2. counts_end: the counts to the end transcript
    3. counts_readthrough: the counts to the readthrough transcript
    4. tail_length_base: a comma-separated list of lengths that the
       readthrough transcript reads extend past the 3' end of the feature
    5. tail_length_end: the tail lengths of the end transcript reads
    6. tail_length_readthrough: the tail lengths of the readthrough transcript
       reads

    SAM file reading follows "Sequence Alignment/Map Format Specification"
    updated April 30, 2020 from the SAM/BAM Format Specification Working Group.

    In addition, if pathToHits is specified, all matching reads are written
    to this path and file. Alignments are written in SAM format (no header
    lines) with 3 added optional fields:
    1. fe:Z:FEATURE_NAME
    2. ba:i:TAIL_LENGTH_BASE
    3. en:i:TAIL_LENGTH_END (or -1 if not an end match)
    4. rt:i:TAIL_LENGTH_READTHROUGH (or -1 if not a readthrough match)
    5. cw:i:1/COUNTING_WEIGHT
    """

    ff=ht.GFF_Reader(pathToGTF,end_included=True)
    exons=ht.GenomicArrayOfSets("auto",stranded=False)
    strands={}
    for feature in ff:
        if feature.type=="exon":
            exons[feature.iv]+=feature.name
            if feature.name not in strands.keys():
                strands[feature.name]=









    #import itertools
    #for feature in itertools.islice(ff,20):
    #    print(feature)

