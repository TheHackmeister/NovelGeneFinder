import sys
import os
import linecache
from itertools import izip, izip_longest
from Hmm import Hmm
from HammerTable import HammerTable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from openFile import openFile

import numpy

# This file should be thought of as the starting point for the program.
# It takes a hammer file, and 2 (paired) fastq files, and splits them into 
# orphans (neither pair matched an HMM) and oeas (one end achored) reads.
# 
# It first loads the hammer table, then opens the fastq files and goes through
# line by line deciding in which file the pair should get put. 
#
# To cut down on overhead, I use write(line, file) to write to files. This 
# method remembers all the files that have been opened, and just reuses those
# pointers. Then when the program is done, closeWrite() closes those files. 
def split(hammerFile, fastqFFile, fastqRFile, flags = {}):
    fOrphanFile = fastqFFile + ".orphans";
    rOrphanFile = fastqRFile + ".orphans";
    fOeaFile = fastqFFile + ".oeas";
    rOeaFile = fastqRFile + ".oeas";

    hmms = HammerTable(hammerFile)

    if flags.get("verbos", False):
        print "Finished hammer table load."
    

    array = [] # For figuring out what the insert size is.
    with openFile(fastqFFile) as fastqf, openFile(fastqRFile) as fastqr:
        for fLine, rLine in izip_longest(fastqf, fastqr):
            if(fLine.startswith("@")):
                fSequence = next(fastqf); 
                rSequence = next(fastqr); 
                if flags["sequenceNames"]: 
                    # If the hammer file has sequences by name rather than 
                    # by genetic code, the format is: 
                    # @TheSequenceName OtherThings
                    # So it basically just removes @ and keeps the sequence 
                    # name.
                    if flags.get("verbos", False):
                        print "fLine:", fLine
                        print "fSeq:", fLine[1:].split()[0]
                    fSeq = hmms.getSequence(fLine[1:].split()[0]) 
                    rSeq = hmms.getSequence(rLine[1:].split()[0])
                else:
                    # Remove the trailing newline.
                    fSeq = hmms.getSequence(fSequence[:-1])
                    rSeq = hmms.getSequence(rSequence[:-1])

                if (not fSeq) and (not rSeq):
                    # Each sequence ends up having 4 lines in the fastq, this 
                    # puts them all in the desired file.
                    write(fLine + fSequence + next(fastqf) + next(fastqf), fOrphanFile)
                    write(rLine + rSequence + next(fastqr) + next(fastqr), rOrphanFile)
                elif fSeq and rSeq:
                    # If the sequences are in DIFFERENT models.
                    if fSeq.getHmm() != rSeq.getHmm():
                        pass # Do something here!
                    else: # Sequences are in the same model.  
                        # To find the insert distance.
                        array.append(rSeq.hmmFrom - fSeq.hmmFrom) 

                else:
                    # Each sequence ends up having 4 lines in the fastq, this 
                    # puts them all in the desired file.
                    write(fLine + fSequence + next(fastqf) + next(fastqf), fOeaFile)
                    write(rLine + rSequence + next(fastqr) + next(fastqr), rOeaFile)

    if flags.get("verbos", False):
        print "Finished Fastq split"


    # This section is mainly used for learning about the data. In production, it
    # will likely not be used.
    if flags.get("createGraphs", False):
        hmms.createGraphs()    
    
    print "Average", numpy.mean(array)
    print "Std", numpy.std(array)
    pyplot.hist(array, bins=numpy.arange(-200,200));
    pyplot.draw();
    pyplot.savefig("insertSizeDistribution.png");

# Writes the line to the file. 
# Because of frequent appending to files, it leaves the file open to save
# overhead. When we're done, closeWrite() should get called to close these.
def write(line, fileName):
    try:
        files[fileName].write(line)
    except KeyError:
        if os.path.exists(fileName):
            raise Exception("The " + fileName + " file already exists.")
        files[fileName] = open(fileName, "a") 
        files[fileName].write(line)

# Closes the files opened by write(line, fileName).
# Should be called at the end of execution, regardless of success.
def closeWrite():
    for f in files:
        files[f].close()

# I wrap in a try block because I'm opening files, but not in a with block. 
# So if anything happens, I can still close those files.
def main():
    try:
        split(sys.argv[1], sys.argv[2], sys.argv[3], {"verbos": True, "sequenceNames": True, "createGraphs": False})
    finally:
        closeWrite()

if __name__ == "__main__":
    files = {} # Global to the file. Used in write(line, fileName) and closeWrite().
    main()


print "Done";
