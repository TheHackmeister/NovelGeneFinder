import sys
from HmmSequence import HmmSequence 
from Hmm import Hmm
from openFile import openFile

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot

# HammerTable represents a hamme table file, along with a couple of helper 
# commands for accessing and displaying the data.
class HammerTable(object):
    def __init__(self, hammerFile, createGraphs = False, calcInsertSize = True):
        self.sequences = {}
        self.hmms = {}
        print "Started processing sequences."
        with openFile(hammerFile) as sequences: # Load the hammer table sequences.
            for line in sequences:
                if line[0] != "#": # Skip commented lines.
                    # l[0] is the sequence.
                    # l[2] is the HMM name.
                    # l[11] is the direction of the sequence.
                    l = line.strip().split()

                    try: 
                        self.sequences[l[0] + l[11]].addHmm(l)
                        # I had added a reference to the hammer table, but I don't think it's needed.
                        # self.sequences[l[0] + l[11]].hammerTable = self
                    except KeyError:
                        self.sequences[l[0] + l[11]] = HmmSequence(line)

                    if createGraphs:
                        try:
                            self.hmms[l[2]].addLine(l) 
                        except KeyError:
                            self.hmms[l[2]] = Hmm(l[2])
                            self.hmms[l[2]].addLine(l)
        print "Finished processing sequences."

        if createGraphs:
            self.createGraphs()

        if calcInsertSize:
            self.calcInsertSize()

    def getHmm(self, hmmName):
        return self.hmms.get(hmmName, False)

    def getSequence(self, sequence):
        print "Sequence:", sequence
        return self.sequences.get(sequence, None)

    def getSequenceHmms(self, sequence):
        try:
            return self.sequences[sequence].getHmms()
        except KeyError:
            return set([])
            
    def createGraphs(self):
        print "Started creating graphs."
        for hmm in self.hmms:
            self.hmms[hmm].drawVectorAndBreakpoints()
            print "Created", hmm
        print "Finished creating graphs."

    def getMissingBlock(self, hmm, forwardName, reverseName):
        return None

    def calcInsertSize(self):
        insertSizeArray = []
        for fSeq in self.sequences:
            if fSeq[:-1] == "-":
                continue # Only need to look at + strands.

            rSeq = self.sequences.get(fSeq[:-1] + "-", None)
            fSeq = self.sequences.get(fSeq, None)
            if not rSeq:
                continue # Only look at strands with a pair.

            # Take the set from the forward and reverse sequences and loops over their intersection.
            for hmm in fSeq.getHmms() & rSeq.getHmms():
                insertSizeArray.append(int(fSeq.getHmm(hmm).hmmTo) - int(rSeq.getHmm(hmm).hmmFrom))
                
        print "Average", numpy.mean(insertSizeArray)
        print "Std", numpy.std(insertSizeArray)
        pyplot.hist(insertSizeArray) #, bins=numpy.arange(-200,200));
        pyplot.draw();
        pyplot.savefig("insertSizeDistribution.png");
                        
# Used for testing the class. 
def main():
    createGraphs = False
    calcInsertSize = True
    print "Started processing hammer file."
    hmms = HammerTable(sys.argv[1], createGraphs, calcInsertSize )
    print "Finished proccessing hammer file."

if __name__ == "__main__":
    main()
