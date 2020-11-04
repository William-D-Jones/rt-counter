def _tails2file(dictTails,path):
    """
    Write the readthrough tails dictionary object to a file.
    """

    fileTails=open(path,mode="w",newline="\n")
    for key in dictTails.keys():
        tails=[str(val) for val in dictTails[key]]
        line="".join([" ".join([key," ".join(tails)]),"\n"])
        fileTails.writelines(line)
    fileTails.close()

