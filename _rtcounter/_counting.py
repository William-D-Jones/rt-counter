import HTSeq as ht
import collections as co

# Parser choices
MODE_CNTALL="count_all"
MODE_CNTNON="count_none"
MODE_CNTFRA="count_fractional"
COUNTMODES=[MODE_CNTALL,MODE_CNTNON,MODE_CNTFRA]

# CIGAR data
MOPS=set(["M","=","X"]) # CIGAR operations interpreted as matches
SOPS=set(["N"]) # CIGAR operations interpreted as skips

# Miscellaneous constants
ERROR_CODES={
        "MULTIMAP":1,
        "ORPHAN":2,
        "FAILCHECKS":3,
        "NOFEAT":4,
        "MULTIFEAT":5,
        "INSMATCH":6}

def overlap(listOfGenomicIntervals):
    """
    Computes length of overlap among a collection of HTSeq
    GenomicInterval objects.
    """

    start=[]
    end=[]
    for iv in listOfGenomicIntervals:
        start.append(iv.start)
        end.append(iv.end)
    ov=min(end)-max(start)

    return ov
        
def _calc_rt(pathToGTF,pathToBAM,pathToHits,pathToFails,fAnchor,rtAnchor,\
        countMultimappers):
    """
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
           overlap any interval lying immediately upstream of the most 
           upstream exon of F.
    2. A template also counts toward the 'readthrough' transcript of F if
       it counts toward F's 'base' transcript and also:
        d. At least one segment in the template has at least 1 CIGAR
           match operations that overlaps an interval lying downstream of the 
           most downstream exon of F.

        For every template satisfying criteria a-d, we calculate the
        'tail length' for each segment:
           The tail length is the length of the interval bounded
           by, but not containing, the most downstream nucleotide of the most 
           downstream exon in F and the most upstream N operation in the
           segment, with respect to F, that is downstream of the most
           downstream match operation in the segment that overlaps F. In
           other words, the interval used for the tail length must contain
           no N operations.

    Returns:

    (i)
    A Pandas dataframe containing these columns, left-to-right:
    0. id: the id of the feature
    1. count_base: the count to the base transcript
    2. count_readthrough: the count to the readthrough transcript

    (ii)
    A dictionary in which each id from the dataframe is a key whose value
    is a list of the tail lengths 

    SAM file reading follows "Sequence Alignment/Map Format Specification"
    updated April 30, 2020 from the SAM/BAM Format Specification Working Group.

    In addition, if pathToHits is specified, all matching reads are written
    to this path and file. Alignments are written in SAM format with these 
    added optional fields:
    1. fe:Z:FEATURE_NAME
    2. ba:i:LEN (count of match operations supporting base transcript)
    3. tl:i:TAIL_LENGTH (or -1 if not a readthrough match)

    Finally, if pathToFails is specified, all nonmatching reads are written
    to this path and file. Alignments are written in SAM format with 1 added 
    optional field:
    1. xx:i:ERROR_CODE
    """

    ff=ht.GFF_Reader(pathToGTF,end_included=True)
    exons=ht.GenomicArrayOfSets("auto",stranded=False)
    ivsExonsDown={} # Holds the most downstream exon of each feature.name
    ivsExonsUp={} # Holds the most upstream exon of each feature.name
    for feature in ff:
        if feature.type=="exon":
            exons[feature.iv]+=feature.name
            if feature.name not in strands.keys():
                strands[feature.name]=feature.iv.strand
            if ((feature.name not in exonsTerm.keys()) or \
                    (feature.iv.strand=="-" and \
                    feature.iv.start<exonsTerm[feature.name].iv.start) or \
                    (feature.iv.strand=="+" and \
                    feature.iv.end>exonsTerm[feature.name].iv.end)):
                if feature.iv.strand="+":
                    start=feature.iv.end-fAnchor
                    end=feature.iv.end
                elif feature.iv.strand="-":
                    start=feature.iv.start
                    end=start+fAnchor
                if feature.iv.strand!=".":
                    ivTerm=ht.GenomicInterval(
                            feature.iv.chrom,start,end,feature.iv.strand)
                    if ivTerm.length>=fAnchor:
                        exonsTerm[feature.name]=ivTerm

    cntBase=co.Counter()
    cntEnd=co.Counter()
    cntRt=co.Counter()
    
    aa=ht.SAM_Reader(pathToBAM)
    if pathToHits is not None:
        w_hits=ht.BAM_Writer.from_BAM_Reader(pathToHits,aa)
    if pathToFails is not None:
        w_fails=ht.BAM_Writer.from_BAM_Reader(pathToFails,aa)
    for bundle in ht.pair_SAM_alignments(aa,bundle=True):

        weightMultimappers=0
        failchecks=False

        # Determine if the template is a multimapper, and assign weight
        if len(bundle)>1:
            if countMultimappers==MODE_CNTALL:
                weightMultimappers=1
            elif countMultimappers==MODE_CNTFRA:
                weightMultimappers=1/len(bundle)
            elif countMultimappers==MODE_CNTNON:
                cntBase["_multimapper"]+=1
                for pair in bundle:
                    writeBAMwithOpts(
                            w_fails,pair,[["xx","i",ERROR_CODES["MULTIMAP"]]])
                continue
            else:
                raise ValueError("Unknown choice for countMultimappers mode.")
        else:
            weightMultimappers=1

        for pair in bundle:

            weightOverlappers=0

            # Check that the alignment is paired
            if None in pair:
                cntBase["_orphan_alignment"]+=weightMultimappers
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["ORPHAN"]]])
                continue

            # Check that all segments are aligned and properly paired
            for seg in pair:
                if ((not seg.aligned) or \
                        seg.failed_platform_qc or \
                        seg.pcr_or_optical_duplicate or \
                        (not seg.proper_pair)):
                    failchecks=True
                    break
            if failchecks:
                failchecks=False
                cntBase["_fail_checks"]+=weightMultimappers
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["FAILCHECKS"]]])
                continue

            # Evaluate counts to the base transcript
            # Identify features and match operation counts supporting each one
            listCntFeat=[]
            featSegs=set()
            listIvFlank=[]
            isMultifeat=False
            for seg in pair:
                collect=True # Should we overrwrite the flank
                ivFlank=None # Holds the latest flank
                cntFeat=co.Counter()
                cigar=seg.cigar
                for op in cigar:
                    if op.type in MOPS:
                        for iv,ff in exons[op.ref_iv].steps():
                            for f in ff:
                                cntFeat[f]+=iv.length
                        if collect:
                            ivFlank=op.ref_iv
                            if seg.iv.strand=="-":
                                collect=False
                featSegs |= set(cntFeat.keys())
                if len(featSegs)>1:
                    isMultifeat=True
                    break
                listCntFeat.append(cntFeat)
                listIvFlank.append(ivFlank)
            # Check if the segments match to any feature
            if len(featSegs)==0:
                cntBase["_no_feature"]+=weightMultimappers
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["NOFEAT"]]])
                continue
            # Check if any segment overlaps multiple features
            if isMultifeat:
                cntBase["_multifeature_segments"]+=weightMultimappers
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["MULTIFEAT"]]])
                continue
            # Check if the match length if sufficient to declare a base match
            isLenPass=False
            lenBase=0
            featBase=list(featSegs)[0]
            isAllSegsMatch=True
            for cntFeat in listCntFeat:
                if cntFeat[featBase]>lenBase:
                    lenBase=cntFeat[featBase]
                if cntFeat[featBase]==0:
                    isAllSegsMatch=False
            if lenBase<fAnchor:
                cntBase["_insufficient_match"]+=weightMultimappers
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["INSMATCH"]]])
                continue
            # Construct the inner template interval
            itiStart=None
            itiEnd=None
            itiChrom=None
            for iv in listIvFlank:
                if ((itiStart is None) or iv.start<itiStart):
                    itiStart=iv.start
                if ((itiEnd is None) or iv.end>itiEnd):
                    itiEnd=iv.end
                if itiChrom is None:
                    itiChrom=iv.chrom
            iti=ht.GenomicInterval(itiChrom,itiStart,itiEnd,".")

            # Record the match to the base transcript
            cntBase[featBase]+=weightMultimappers

            # Evaluate counts to the end transcript
            listLenEndMatch=[]
            if not nucsTerm[featBase].overlaps(iti):
                ivEndFeat=None
                for seg in pair:
                    listLenEndMatch.append(0)
                    lenEndMatch=0
                    cigar=seg.cigar
                    for op in cigar:
                        if op.type in MOPS:
                            lenEndMatch+=op.size
                            if ivEndFeat is None:
                                ivEndFeat=op.ref_iv
                            else:
                                ivEndFeat.extend_to_include(op.ref_iv)
                        elif op.type in SOPS:
                            if ivEndFeat.overlaps(nucsTerm[featBase]):
                                listLenEndMatch[-1]=lenEndMatch
                                break
                            else:
                                ivEndFeat=None
            else:
                listLenEnd.append(overlap([iit,exonsTerm[featBase]))



            if not isEnd:
                writeBAMwithOpts(w_hits,pair,
                        [["fe","Z",featBase],
                            ["ba","i",lenBase],
                            ["en","i",0]])
                continue
            if isEnd:
                if max(lenEnd)<fAnchor:
                    writeBAMwithOpts(w_hits,pair,
                            [["fe","Z",featBase],
                                ["ba","i",lenBase],
                                ["en","i",0]])
                    continue
                    



                





    if pathToHits is not None:
        w_hits.close()
    if pathToFails is not None:
        w_fails.close()

def writeBAMwithOpts(writer,pair,opts):
    """
    Writes a pair of BAM alignment with the HTSeq BAM writer, 
    adding optional fields.
    writer: an HTSeq BAM writer
    pair: an iterable of HTSeq SAM alignments
    opts: a list of optional fields in the format [[TAG,TYPE,VALUE],...]
        The TAG and TYPE elements must be strings. If the VALUE element is
        not a string already, it will be converted to a string.
    """

    # Convert optional field VALUE to string if needed
    for i in range(len(opts)):
        opts[i][2]=str(opts[i][2])

    for aln in pair:

        # Construct the SAM alignment with new optional fields
        line=aln.get_sam_line()
        for opt in opts:
            line="\t".join([line,":".join(opt)])
        aln_new=ht.SAM_Alignment.from_SAM_line(line)

        # Write the line
        writer.write(aln_new)
