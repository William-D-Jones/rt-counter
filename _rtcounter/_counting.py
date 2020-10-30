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
ERROR_CODES={"ORPHAN":"1","FAILCHECKS":"2","INSMATCH":3,"NOFEAT":4,
        "MULTIMAP":5,"MULTIFEAT":6,"MULTIFEAT_INT":7}

def _calc_rt(pathToGTF,pathToBAM,pathToHits,pathToFails,fAnchor,rtAnchor,\
        countMultimappers):
    """
    Counts paired-end reads from a BAM file using the following procedure:
    1. A template counts toward a feature's 'base' transcript if:
        a. At least one segment in the template has at least fAnchor CIGAR
           match operations (M, =, or X) that match a feature F.
        b. No nucleotides in either segment of the template match any features
           other than F.
        c. The 'inner template interval' contains no features other than F.
           The inner template interval is the smallest genomic interval
           containing (in other words, flanked by):
            i.  The leftmost M operation of the + strand segment for which
                no N operations lie between this M operation and the
                rightmost CIGAR operation in the segment.
            ii. The rightmost M operation of the - strand segment for 
                which no N operations lie between this M operation and 
                the leftmost CIGAR operation in the segment.
    2. A template counts toward a feature's 'end' transcript if:
        a. It counts toward the feature's 'base' transcript.
        b. One of the following is true:
            i.  The inner template interval contains at least fAnchor 
                nucleotides of the most 3' exon in F including the most 3'
                nucleotide of this exon, or
            ii. Condition (i) is false but at least one of the segments
                contains at least fAnchor match operations, uninterrupted by
                any N operations, that match the most 3' exon in F and the
                most 3' nucleotide of this exon.
    3. A template counts toward a feature's 'readthrough' transcript if:
        a. It counts toward the feature's 'end' transcript.
        b. The 'tail length' is at least rtAnchor. For every template
           matching the end transcript, the tail length is the distance of 
           the path starting from the most 3' nucleotide of the most 3' exon 
           of the feature to the end of the inner template interval such
           that this path follows the same direction of transcription as the
           feature.

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
    4. tail_length_end: the tail lengths of the end transcript reads
    5. tail_length_readthrough: the tail lengths of the readthrough transcript
       reads

    SAM file reading follows "Sequence Alignment/Map Format Specification"
    updated April 30, 2020 from the SAM/BAM Format Specification Working Group.

    In addition, if pathToHits is specified, all matching reads are written
    to this path and file. Alignments are written in SAM format with 3 added 
    optional fields:
    1. fe:Z:FEATURE_NAME
    2. ba:i:LEN (count of match operations supporting base transcript)
    3. en:i:LEN (count of match operations supporting end transcript)
    4. te:i:TAIL_LENGTH_END (-1 if not an end match)
    5. rt:i:LEN (count of match operations supporting readthrough transcript)
    6. tr:i:TAIL_LENGTH_READTHROUGH (-1 if not a readthrough match)
    7. cw:i:1/COUNTING_WEIGHT

    Finally, if pathToFails is specified, all nonmatching reads are written
    to this path and file. Alignments are written in SAM format with 1 added 
    optional field:
    1. xx:i:ERROR_CODE
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
                if len(cntFeat.keys())>1 or len(featSegs)>1:
                    isMultifeat=True
                    break
                listCntFeat.append(cntFeat)
                listIvFlank.append(ivFlank)
            # Check if any segment overlaps multiple features
            if isMultifeat:
                cntBase["_multifeature_segments"]+=weightMultimappers
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["MULTIFEAT"]]])
                continue
            # Construct and check the inner template interval
            isItiMultifeat=False
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
            itiFeat=set()
            for iv,ff in exons[iti].steps():
                itiFeat |= ff
            if len(itiFeat)>1:
                isItiMultifeat=True
            # Check if the iti overlaps multiple features
            if isItiMultifeat:
                cntBase["_multifeature_template"]+=weightMultimappers
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["MULTIFEAT_INT"]]])
                continue
            # Check if the match length if sufficient to declare a base match
            isLenPass=False
            lenBase=0
            if len(featSegs)==0:
                cntBase["_no_feature"]+=weightMultimappers
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["NOFEAT"]]])
                continue
            featBase=list(featSegs)[0]
            for cntFeat in listCntFeat:
                if cntFeat[featBase]>lenBase:
                    lenBase=cntFeat[featBase]
            if lenBase<fAnchor:
                cntBase["_insufficient_match"]+=weightMultimappers
                writeBAMwithOpts(
                        w_fails,pair,[["xx","i",ERROR_CODES["INSMATCH"]]])
                continue

            # Record the match to the base transcript
            writeBAMwithOpts(w_hits,pair,
                    [["fe","Z",featBase],["ba","i",lenBase]])
            cntBase[featBase]+=1

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
