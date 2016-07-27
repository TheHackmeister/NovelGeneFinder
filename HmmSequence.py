from openFile import openFile

# This represents a genetic sequence from an HMM model.
class HmmSequence(object):
    def __init__(self, line):
        l = line.split();
        self.sequence = l[0];
#        print "Sequence:", self.sequence
        self.hmmName = l[2];
        self.hmmFrom = int(l[4]);
        self.hmmTo = int(l[5]);
        # Commenting these out to save space.
#        self.alignFrom = int(l[6]);
#        self.alignTo = int(l[7]);
        self.sequenceLength = int(l[10]);
 #       self.strand = l[11];

        self.alignedLength = abs(int(l[7]) - int(l[6])) # l[6/7] are alignFrom and alignTo. 

    def __str__(self):
        return self.hmmName;

    
    def isSplit(self, percent = .8):
        return float(self.alignedLength)/float(self.sequenceLength) <= percent

    # Tries to get the sequence's hmm.
    # Returns false if it doesn't exist.
    def getHmm(self):
        try:
            return self.hammerTable.getHmm(self.hmmName)
        except NameError: 
            return False

    # I don't think this is used. Will remove eventually.
    @staticmethod
    def createHmmSequences(fileLoc):
        sequenceList = [];
        with openFile(fileLoc) as sequences:
            for line in sequences:
                if line[0] != "#":
                    sequenceList.append(HmmSequence(line));
        return sequenceList;
