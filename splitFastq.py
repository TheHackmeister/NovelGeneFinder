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


    basePath = os.path.dirname(fastqFFile)

    fOrphanFile = fastqFFile + ".orphans"
    rOrphanFile = fastqRFile + ".orphans"
    fOeaFile = fastqFFile + ".oeas"
    rOeaFile = fastqRFile + ".oeas"

    hmms = HammerTable(hammerFile, calcInsertSize = False)

    
    if flags.get("verbos", False):
        print "Finished hammer table load."
    

    with openFile(fastqFFile) as fastqf, openFile(fastqRFile) as fastqr:
        for fLine, rLine in izip_longest(fastqf, fastqr):
            if(fLine.startswith("@")):
                fSequence = next(fastqf) 
                rSequence = next(fastqr) 
                if flags["sequenceNames"]: 
                    # If the hammer file has sequences by name rather than 
                    # by genetic code, the format is: 
                    # @TheSequenceName OtherThings
                    # So it basically just removes @, keeps the sequence 
                    # name, and adds a strand direction (+/-).
                    forwardName = fLine[1:].split()[0] + "+"
                    reverseName = rLine[1:].split()[0] + "-"
                else:
                    # Remove the trailing newline.
                    forwardName = fSequence[:-1]
                    reverseName = rSequence[:-1]
                
                fHmms = hmms.getSequenceHmms(forwardName) 
                rHmms = hmms.getSequenceHmms(reverseName)

                combinedHmms = fHmms | rHmms # The union of the fHmms and rHmms

                # These are orphan reads.
                if len(combinedHmms) == 0:
                    # Each sequence ends up having 4 lines in the fastq, this 
                    # puts them all in the desired file.
                    #write(fLine + fSequence + next(fastqf) + next(fastqf), fOrphanFile)
                    #write(rLine + rSequence + next(fastqr) + next(fastqr), rOrphanFile)
                    pass
                else:
                    for hmm in combinedHmms:
                        # Create one big file with all the sequences. 
                        if flags.get("completeSequenceFile", False):
                            write(fLine + fSequence + next(fastqf) + next(fastqf), os.path.join(basePath, hmm, hmm + "-full-forward.fastq"))
                            write(rLine + rSequence + next(fastqr) + next(fastqr), os.path.join(basePath, hmm, hmm + "-full-backward.fastq"))

                        # They both map to the HMM. Don't need to do anything. 
                        if hmm in fHmms and hmm in rHmms:
                            pass

                        # Only the forward strand maps to the HMM. 
                        elif hmm in fHmms:
                            block = hmms.getMissingBlock(hmm, forwardName, reverseName)
                            if block != None:
                                write(rLine + rSequence + next(fastqr) + next(fastqr), os.path.join(basePath, hmm, hmm + "-block-" + block + ".fastq"))
                                next(fastqf) # I moved forward the reverse file, so I have to do the same to the forward file.
                                next(fastqf)

                        # Only the reverse strand maps to the HMM.
                        else:
                            block = hmms.getMissingBlock(hmm, forwardName, reverseName)
                            if block != None:
                                write(fLine + fSequence + next(fastqf) + next(fastqf), os.path.join(basePath, hmm, hmm + "-block-" + block + ".fastq"))
                                next(fastqr) # I moved forward the forward file, so I have to do the same to the reverse file.
                                next(fastqr)

    if flags.get("verbos", False):
        print "Finished Fastq split"


# Writes the line to the file. 
# Because of frequent appending to files, it leaves the file open to save
# overhead. When we're done, closeWrite() should get called to close these.
def write(line, fileName):
    if not os.path.exists(os.path.dirname(fileName)):
        try:
            os.makedirs(os.path.dirname(fileName))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
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
        split(sys.argv[1], sys.argv[2], sys.argv[3], {"verbos": True, "sequenceNames": True, "createGraphs": False, "completeSequenceFile": True})
    finally:
        closeWrite()

if __name__ == "__main__":
    files = {} # Global to the file. Used in write(line, fileName) and closeWrite().
    main()


print "Done";
