import sys
from HmmSequence import HmmSequence 
from Hmm import Hmm
from openFile import openFile

# HammerTable represents a hamme table file, along with a couple of helper 
# commands for accessing and displaying the data.
class HammerTable(object):
    def __init__(self, hammerFile):
        self.sequences = {}
        with openFile(hammerFile) as sequences: # Load the hammer table sequences.
                for line in sequences:
                    if line[0] != "#": # Skip commented lines.
                        newSequence = HmmSequence(line)
                        newSequence.hammerTable = self
                        self.sequences[newSequence.sequence] = newSequence
        self.hmms = Hmm.createHmmsFromHammer(hammerFile) # Create the hammer table HMM models.
        self.hmms = {}

    def getHmm(self, hmmName):
        return self.hmms.get(hmmName, False)

    def getSequence(self, sequence):
        return self.sequences.get(sequence, None)
            
    def createGraphs(self):
        for hmm in self.hmms:
            self.hmms[hmm].drawVectorAndBreakpoints()
            print "Created", hmm
                        
# Used for testing the class. 
def main():
    HammerTable(sys.argv[1])

if __name__ == "__main__":
    main()
