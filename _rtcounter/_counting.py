import HTSeq as ht
import collections as co
import datetime as dt

# CIGAR data
MOPS=set(["M","=","X"]) # CIGAR operations interpreted as matches
SOPS=set(["N"]) # CIGAR operations interpreted as skips
ROPS=set(["M","=","X","N","D"]) # Reference-consuming CIGAR operations

# Miscellaneous constants
ERROR_CODES={
        "MULTI_MAP":1,
        "ORPHAN":2,
        "FAIL_CHECKS":3,
        "NO_NAME":4,
        "MULTI_NAME":5,
        "INSUFF_MATCH":6,
        "UPSTREAM_MATCH":7}

def msg(message):
    """
    Print a string with date and time.
    """
    txt="\t".join([dt.datetime.now(),message])
    print(txt)

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
        
def _calc_rt(pathToGTF,pathToBAM,pathToHits,pathToFails,fAnchor):
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
    is a list of the tail lengths.

    SAM file reading follows "Sequence Alignment/Map Format Specification"
    updated April 30, 2020 from the SAM/BAM Format Specification Working Group.

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
    """

    msg("Reading feature file...")
    ff=ht.GFF_Reader(pathToGTF,end_included=True)
    exons=ht.GenomicArrayOfSets("auto",stranded=False)
    # The ivBounds dictionary has feature names as its keys and values
    # consisting of genomic intervals encompassing all exons with each name,
    # but only exons that have an associated strand + or -
    ivBounds={}
    for feature in ff:
        if feature.type=="exon":
            exons[feature.iv]+=feature.name
            if feature.iv.strand!=".":
                if (feature.name not in ivBounds.keys()):
                    ivBounds[feature.name]=feature.iv
                else:
                    ivBounds[feature.name].extend_to_include(feature.iv)
    msg("".join(["Read ",str(len(exons.keys()))," named groups of exons."]))

    cntBase=co.Counter() # Counter of base transcripts and errors
    cntRT=co.Counter() # Counter of readthrough transcripts
    dictRT=co.defaultdict(list) # Lists of tail lengths for each feature name
    prog={"READS":0} # Holds progress statistics
    
    aa=ht.SAM_Reader(pathToBAM)
    if pathToHits is not None:
        w_hits=ht.BAM_Writer.from_BAM_Reader(pathToHits,aa)
    else:
        w_hits=None
    if pathToFails is not None:
        w_fails=ht.BAM_Writer.from_BAM_Reader(pathToFails,aa)
    else:
        w_fails=None

    msg("Started counting reads...")

    for bundle in ht.pair_SAM_alignments(aa,bundle=True):

        # Progress update
        prog["READS"]+=1
        if prog["READS"]>0 and prog["READS"]%10000==0:
            msg("".join(["In progress: Counted ",str(prog["READS"]),
                " reads, ",str(cntBase["_HITS"])," base hits, and ",
                str(cntRT["_HITS"])," readthrough hits."]))

        # Reject multimappers
        if len(bundle)>1:
            cntBase["_MULTI_MAP"]+=1
            for pair in bundle:
                writeBAMwithOpts(
                    w_fails,pair,[["xx","i",ERROR_CODES["MULTI_MAP"]]])
            continue

        for pair in bundle:

            # Check that the alignment is paired
            if None in pair:
                cntBase["_ORPHAN"]+=1
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["ORPHAN"]]])
                continue

            # Check that all segments are aligned and properly paired
            failchecks=False
            for seg in pair:
                if ((not seg.aligned) or \
                        seg.failed_platform_qc or \
                        seg.pcr_or_optical_duplicate or \
                        (not seg.proper_pair)):
                    failchecks=True
                    break
            if failchecks:
                failchecks=False
                cntBase["_FAIL_CHECKS"]+=1
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["FAIL_CHECKS"]]])
                continue

            # Identify a unique feature name for the template, if it exists
            # Identify feature names and the supporting match operation counts
            listCntNames=[]
            namesBase=set()
            isMultiName=False # Set to True if an op overlaps multiple names
            for seg in pair:
                cntNames=co.Counter()
                for op in seg.cigar:
                    if op.type in MOPS:
                        for iv,names in exons[op.ref_iv].steps():
                            if len(names)>1:
                                isMultiName=True
                                break
                            for name in names:
                                cntNames[name]+=iv.length
                    if isMultiName:
                        break
                if isMultiName:
                    break
                namesBase |= set(cntNames.keys())
                listCntNames.append(cntNames)
            # Check if any segment overlaps multiple feature names
            if isMultiName or len(namesBase)>1:
                cntBase["_MULTI_NAME"]+=1
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["MULTI_NAME"]]])
                continue
            # Check if any segment overlaps any feature name
            if len(namesBase)==0:
                cntBase["_NO_NAME"]+=1
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["NO_NAME"]]])
                continue
            # Check if the overlap length is sufficient to declare a base match
            nameBase=list(namesBase)[0]
            lenBase=max(list(listCntNames.values()))
        if lenBase<fAnchor:
            cntBase["_INSUFF_MATCH"]+=1
            writeBAMwithOpts(
                    w_fails,pair,[["xx","i",ERROR_CODES["INSUFF_MATCH"]]])
            continue

        # Extract information about the matched name
        ivBase=ivBounds[nameBase]
        sBase=ivBase.strand

        # Identify upstream matches (upstream of the most upstream exon)
        isUp=False # Set to True if we find an upstream match
        for seg in pair:
            if ((sBase=="+" and seg.iv.start<ivBase.start) or \
                    (sBase=="-" and seg.iv.end>ivBase.end)):
                isUp=True
                break
        if isUp:
            cntBase["_UPSTREAM_MATCH"]+=1
            writeBAMwithOpts(
                    w_fails,pair,
                    [["xx","i",ERROR_CODES["UPSTREAM_MATCH"]]])
            continue

        # Identify readthrough (matches downstream of the most downstream exon)
        for seg in pair:

            # Check if the segment spans downstream of the exons at all
            if ((sBase=="+" and seg.iv.end<=ivBase.end) or \
                    (sBase=="-" and seg.iv.start>=ivBase.start)):
                continue
            
            # Construct an interval of length 0, to be included to include
            # any CIGAR operations that qualify as readthrough
            if sBase=="+":
                posBound=ivBase.end
            elif sBase=="-":
                posBound=ivBase.start
            ivRT=ht.GenomicInterval(
                    ivBase.chrom,posBound,posBound,sBase)

            # Construct a downstream iterator through the CIGAR
            if seg.iv.strand=="+":
                itCIGAR=seg.cigar
            elif seg.iv.strand=="-":
                itCIGAR=reversed(seg.cigar)

            # Update ivRT to include any readthrough
            isReading=False # True if we are in readthrough operations
            for op in itCIGAR:
                if op in ROPS:
                    if op not in SOPS:
                        isReading=True
                        if ((sBase=="+" and op.ref_iv.end>ivBase.end) or \
                                (sBase=="-" and op.ref_iv.start<ivBase.start)):
                            ivRT.extend_to_include(op.ref_iv)
                    elif isReading: # Encountered a skip, so stop counting

                        # Correct the upstream bound of ivRT
                        if sBase=="+":
                            ivRT.start=ivBase.end
                        elif sBase=="-":
                            ivRT.end=ivBase.start

                        listLenRT.append(ivRT.length)
                        break
        lenRT=max(listLenRT)

        # Increment counters
        cntBase[nameBase]+=1
        cntBase["_HITS"]+=1
        if lenRT>0:
            cntRT[nameBase]+=1
            cntRT["_HITS"]+=1
            dictRT[nameBase].append(lenRT)

        # Write the hit
        if pathToHits is not None:
            opts=[["fe","Z",nameBase],["ba","i",lenBase]]
            if lenRT>0:
                opts.append(["tl","i",lenRT])
            writeBAMwithOpts(w_hits,pair,opts)

    if pathToHits is not None:
        w_hits.close()
    if pathToFails is not None:
        w_fails.close()

    msg("".join["Finished counting reads. Summary of results:",
        "Total reads from file: ",str(prog["READS"]),
        "Total base hits: ",str(cntBase["_HITS"]),
        "Total readthrough hits: ",str(cntRT["_HITS"]),
        "Total multi-mapping reads: ",str(cntBase["_MULTI_MAP"]),
        "Total orphan reads: ",str(cntBase["_ORPHAN"]),
        "Total reads failing checks: ",str(cntBase["_FAIL_CHECKS"]),
        "Total reads matching no feature name: ",str(cntBase["_NO_NAME"]),
        "Total reads matching multiple names: ",str(cntBase["_MULTI_NAME"]),
        "Total reads with too few matches: ",str(cntBase["_INSUFF_MATCH"]),
        "Total reads with upstream match: ",str(cntBase["_UPSTREAM_MATCH"])])

    # Need to return objects

def writeBAMwithOpts(writer,pair,opts):
    """
    Writes a pair of BAM alignment with the HTSeq BAM writer, 
    adding optional fields.
    writer: an HTSeq BAM writer. If None, we do nothing.
    pair: an iterable of HTSeq SAM alignments
    opts: a list of optional fields in the format [[TAG,TYPE,VALUE],...]
        The TAG and TYPE elements must be strings. If the VALUE element is
        not a string already, it will be converted to a string.
    """

    if writer is None:
        return

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
