import HTSeq as ht
import collections as co

# Parser choices
MODE_CNTALL="count_all"
MODE_CNTNON="count_none"
MODE_CNTFRA="count_fractional"
COUNTMODES=[MODE_CNTALL,MODE_CNTNON,MODE_CNTFRA]

# CIGAR data
MOPS=set("M","=","X") # CIGAR operations interpreted as matches
SOPS=set("N") # CIGAR operations interpreted as skips

def _calc_rt(pathToGTF,pathToBAM,pathToHits,pathToFails,fAnchor,rtAnchor,\
        countMultimappers,countOverlappers):
    """
    Counts paired-end reads from a BAM file using the following procedure:
    1. A template counts toward a feature's 'base' transcript if:
        a. At least one segment in the template has at least fAnchor CIGAR
           match operations on its 3' end, uninterrupted by any skip (N)
           operations, that matches a feature F.
        b. No nucleotides in either segment of the template match any features
           other than F.
        c. The other segment in the template has at least fAnchor CIGAR match
           operations on its 3' end uninterrupted by any skip (N) operations,
           that may or may not match F.
        d. The genomic interval between the most 3' non-N operation of the
           first segment and the most 3' non-N operation of the segment segment
           contains no features other than F.
    2. A template counts toward a feature's 'end' transcript if:
        a. It counts toward the feature's 'base' transcript.
        b. The genomic interval specified in (1d) contains fAnchor of the most
           3' nucleotides of the most 3' exon in F.
    3. A template counts toward a feature's 'readthrough' transcript if:
        a. It counts toward the feature's 'end' transcript.
        b. The genomic interval specified in (1d) contains rtAnchor of the
           nucleotides immediately following the most 3' exon in F.

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
    exonsLatest={} # Holds the most 3' exon of each feature.name so far
    strands={}
    for feature in ff:
        if feature.type=="exon":
            exons[feature.iv]+=feature.name
            if feature.name not in strands.keys():
                strands[feature.name]=feature.iv.strand
            if ((feature.name not in exonsLatest.keys()) or \
                    (feature.iv.strand=="-" and \
                    feature.iv.end<exonsLatest[feature.name].iv.end) or \
                    (feature.iv.strand=="+" and \
                    feature.iv.end>exonsLatest[feature.name].iv.end)):
                exonsLatest[feature.name]=feature
    exonsTerm=ht.GenomicArrayOfSets("auto",stranded=False)
    for key in exonsLatest.keys():
        feature=exonsLatest[key]
        exonsTerm[feature.iv]+=feature.name

    cntBase=co.Counter()
    cntEnd=co.Counter()
    cntRt=co.Counter()
    
    aa=ht.SAM_Reader(pathToBAM)
    for alns in ht.pair_SAM_alignments(aa,bundle=True):

        weightMultimappers=0
        unmapped=False

        # Determine if the template is a multimapper, and assign weight
        if len(alns)>1:
            if countMultimappers==MODE_CNTALL:
                weightMultimappers=1
            elif countMultimappers==MODE_CNTFRA:
                weightMultimappers=1/len(bundle)
            elif countMultimappers==MODE_CNTNON:
                cntBase["_multimapper"]+=1
                continue
            else:
                raise ValueError("Unknown choice for countMultimappers mode.")
        else:
            weightMultimappers=1

        for segs in alns:

            weightOverlappers=0

            # Check that the alignment is paired
            if None in segs:
                cntBase["_orphan_alignment"]+=weightMultimappers
                continue

            # Check that all segments are mapped
            for seg in segs:
                if not seg.aligned:
                    unmapped=True
                    break
            if unmapped:
                unmapped=False
                cntBase["_unmapped"]+=weightMultimappers
                continue

            # Identify base transcripts to which the templates matches

            # First, identify candidate matches

            # The 'match region' is the portion of the CIGAR that we require to 
            # contain a match, that is, the longest possible continuous portion 
            # of the CIGAR containing the CIGAR's rightmost match operation 
            # (M, =, or X) but not containing any skips (N)

            id_all=[] # List of ID sets matching each segment
            id_match=[] # List of ID sets matching each segment in match region
            len_match_any=[] # Length of match region in each segment
            len_match_feat=[] # Length of match region matching a feature
            for i in range(len(segs)):

                len_match.append(0)
                id_all.append(set())
                id_match.append(set())

                cigar=seg.cigar
                inMatchReg=True # Tracks if we have encountered an N

                for j in reversed(range(len(cigar))):
                    if cigar[j].type not in MOPS:
                        if cigar[j].type in SOPS:
                            inMatchReg=False
                        continue
                    for iv,ff in exons[cigar[j].ref_iv].steps():
                        id_all[i] |= ff
                        if inMatchReg:
                            id_match[i] |= ff
                            len_match_any[i]+=cigar[j].size
                            if len(ff)>0:
                                len_match_feat[i]+=cigar[j].size

            # Second, decide if any of the candidate matches constitutes a
            # match to the corresponding 'base' transcript

            ids=set()
            isLenPassAny=True
            isLenPassFeat=False
            for i in len(id_all):
                ids |= id_all[i]
                if len_match_feat[i]>=fAnchor:
                    isLenPassFeat=True
                if len_match_any[i]<fAnchor:
                    isLenPassAny=False
            if not isLenPassFeat:
                continue
            if not isLenPassAny:
                continue
            cntBase[list(id_all)[0]]+=1


            
