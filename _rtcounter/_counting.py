import HTSeq as ht
import collections as co
import datetime as dt
import pandas as pd

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
        "UPSTREAM_MATCH":7,
        "DOWNSTREAM_MATCH":8}

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

    _msg("Reading feature file...")
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
    _msg("".join(["Read ",str(len(ivBounds.keys()))," named sets of exons."]))

    cntStats=co.Counter() # Counter of template statistics
    cntBase=co.Counter() # Counter of base transcripts
    cntRT=co.Counter() # Counter of readthrough transcripts
    dictTails=co.defaultdict(list) # Lists of tail lengths for each feature name
    
    aa=ht.SAM_Reader(pathToBAM)
    if pathToHits is not None:
        w_hits=ht.BAM_Writer.from_BAM_Reader(pathToHits,aa)
    else:
        w_hits=None
    if pathToFails is not None:
        w_fails=ht.BAM_Writer.from_BAM_Reader(pathToFails,aa)
    else:
        w_fails=None

    _msg("Started counting reads...")

    for bundle in ht.pair_SAM_alignments(aa,bundle=True):

        # Progress update
        cntStats["TOT_PAIRS"]+=1
        if cntStats["TOT_PAIRS"]>0 and cntStats["TOT_PAIRS"]%100000==0:
            _msg("".join(["In progress: Counted ",str(cntStats["TOT_PAIRS"]),
                " read pairs, ",str(cntStats["TOT_BASE"])," base hits, and ",
                str(cntStats["TOT_RT"])," readthrough hits."]))

        # Reject multimappers
        if len(bundle)>1:
            cntStats["MULTI_MAP"]+=1
            for pair in bundle:
                writeBAMwithOpts(
                    w_fails,pair,[["xx","i",ERROR_CODES["MULTI_MAP"]]])
            continue

        for pair in bundle:

            # Check that the alignment is paired
            if None in pair:
                cntStats["ORPHAN"]+=1
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
                cntStats["FAIL_CHECKS"]+=1
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
                cntStats["MULTI_NAME"]+=1
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["MULTI_NAME"]]])
                continue
            # Check if any segment overlaps any feature name
            if len(namesBase)==0:
                cntStats["NO_NAME"]+=1
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["NO_NAME"]]])
                continue
            # Check if the overlap length is sufficient to declare a base match
            nameBase=list(namesBase)[0]
            lenBase=0
            for cntNames in listCntNames:
                if cntNames[nameBase]>lenBase:
                    lenBase=cntNames[nameBase]
            if lenBase<fAnchor:
                cntStats["INSUFF_MATCH"]+=1
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
                cntStats["UPSTREAM_MATCH"]+=1
                writeBAMwithOpts(
                        w_fails,pair,
                        [["xx","i",ERROR_CODES["UPSTREAM_MATCH"]]])
                continue

            # Identify readthrough (downstream of the most downstream exon)
            listLenRT=[0]
            isDown=False # Set to True if we find a downstream match
            for seg in pair:

                # Check if the segment spans downstream of the exons at all
                if ((sBase=="+" and seg.iv.end<=ivBase.end) or \
                        (sBase=="-" and seg.iv.start>=ivBase.start)):
                    continue

                isDown=True
                
                # Construct an interval of length 0, to be included to include
                # any CIGAR operations that qualify as readthrough
                if sBase=="+":
                    posBoundOut=ivBase.end
                elif sBase=="-":
                    posBoundOut=ivBase.start
                ivRT=ht.GenomicInterval(
                        ivBase.chrom,posBoundOut,posBoundOut,".")

                # Construct a downstream iterator through the CIGAR
                if sBase=="+":
                    itCIGAR=seg.cigar
                elif sBase=="-":
                    itCIGAR=reversed(seg.cigar)

                # Update ivRT to include any readthrough
                isReading=False # True if we are in readthrough operations
                for op in itCIGAR:
                    if op.type in ROPS:
                        if ((sBase=="+" and op.ref_iv.end>ivBase.end) or \
                                (sBase=="-" and op.ref_iv.start<ivBase.start)):
                            isReading=True
                            if op.type not in SOPS:
                                ivRT.extend_to_include(op.ref_iv)
                        if (isReading and (op.type in SOPS)): # Stop counting
                            # Correct the upstream bound of ivRT
                            if sBase=="+":
                                ivRT.start=ivBase.end
                            elif sBase=="-":
                                ivRT.end=ivBase.start
                            listLenRT.append(ivRT.length)
                            break
            lenRT=max(listLenRT)

            # Check if we have a downstream match but not a readthrough hit
            if isDown and lenRT==0:
                cntStats["DOWNSTREAM_MATCH"]+=1
                writeBAMwithOpts(
                        w_fails,pair,
                        [["xx","i",ERROR_CODES["DOWNSTREAM_MATCH"]]])
                continue

            # Increment counters
            cntBase[nameBase]+=1
            cntStats["TOT_BASE"]+=1
            if lenRT>0:
                cntRT[nameBase]+=1
                cntStats["TOT_RT"]+=1
                dictTails[nameBase].append(lenRT)

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

    _msg("".join(["Finished counting reads.","\n\t","Summary of results:",
        "\n\t","Total reads from file: ",str(cntStats["TOT_PAIRS"]),"\n\t",
        "Total base hits: ",str(cntStats["TOT_BASE"]),"\n\t",
        "Total readthrough hits: ",str(cntStats["TOT_RT"]),"\n\t",
        "Total multi-mapping reads: ",str(cntStats["MULTI_MAP"]),"\n\t",
        "Total orphan reads: ",str(cntStats["ORPHAN"]),"\n\t",
        "Total reads failing checks: ",str(cntStats["FAIL_CHECKS"]),"\n\t",
        "Total reads matching no feature name: ",str(cntStats["NO_NAME"]),
        "\n\t",
        "Total reads matching multiple names: ",str(cntStats["MULTI_NAME"]),
        "\n\t",
        "Total reads with too few matches: ",str(cntStats["INSUFF_MATCH"]),
        "\n\t",
        "Total reads with upstream match: ",str(cntStats["UPSTREAM_MATCH"]),
        "\n\t",
        "Total reads with downstream match but failing readthrough: ",
        str(cntStats["DOWNSTREAM_MATCH"])]))

    # Add zeros to cntRT dictionary to prepare for dataframe
    for key in cntBase.keys():
        if key not in cntRT.keys():
            cntRT[key]=0
    ids={key:key for key in cntBase.keys()}
    # Construct dataframe
    dfCnt=pd.DataFrame(
            {"id":ids,"count_base":cntBase,"count_readthrough":cntRT})
    return dfCnt,dictTails

def _msg(message):
    """
    Print a string with date and time.
    """
    txt="\t".join([str(dt.datetime.now()),message])
    print(txt)

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
