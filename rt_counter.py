import argparse as ap
import csv

# Private modules
import _rtcounter as rt

# Parser defaults
D_FANCHOR=10
D_RTANCHOR=10

# Parse the command line
paMain=ap.ArgumentParser(description=
        "Counts templates showing evidence of transcriptional read-through \
        in paired-end RNA-seq data.")
paMain.add_argument("GTF",type=str,help="Path to Ensembl GTF file containing \
        the genomic features, or any feature file readable by \
        HTSeq.GFF_Reader with end_included=True.",
        metavar="PATH_TO_GTF")
paMain.add_argument("BAM",type=str,help="Path to BAM file containing the \
        aligned reads. BAM alignments must be sorted by read name (QNAME).",
        metavar="PATH_TO_BAM")
paMain.add_argument("pathToCounts",type=str,help="Path and filename to which \
        the output tsv file of counting results will be written.",
        metavar="PATH_TO_OUTPUT_COUNTS")
paMain.add_argument("pathToTails",type=str,help="Path and filename to which \
        the output tsv file of readthrough tail lengths will be written.",
        metavar="PATH_TO_OUTPUT_TAILS")
paMain.add_argument("-m","--pathToHits",type=str,help="If present, writes \
        alignments from the BAM file (no headers) to this path and filename \
        if a match is found to a GTF entry, with 3 added optional fields: \
        fe:Z:FEATURE_NAME, ba:i:tail_length_base, \
        rt:i:tail_length_readthrough",metavar="PATH_TO_HITS")
paMain.add_argument("-x","--pathToFails",type=str,help="If present, writes \
        alignments from the BAM file (no headers) to this path and filename \
        if the alignment is not counted toward any GTF entry, with 1 added \
        option field: xx:i:ERROR_CODE",default=None,metavar="PATH_TO_FAILS")
paMain.add_argument("-f","--featureAnchor",type=int,help="Number of bases in \
        reads that must match the 3' end of each feature in the GTF file. \
        Defaults to "+str(D_FANCHOR)+".",default=D_FANCHOR,metavar=\
        "NUMBER_OF_BASES")
paNames=paMain.parse_args()

# Count read pairs
dfCnt,dictTails=rt._counting._calc_rt(
        pathToGTF=paNames.GTF,
        pathToBAM=paNames.BAM,
        pathToHits=paNames.pathToHits,
        pathToFails=paNames.pathToFails,
        fAnchor=paNames.featureAnchor)

# Write the count dataframe to file
dfCnt.to_csv(paNames.pathToCounts,sep="\t",quoting=csv.QUOTE_NONE,index=False,
        mode="w",line_terminator="\n")

# Write the tails dictionary to file
rt._io._tails2file(dictTails,paNames.pathToTails)

