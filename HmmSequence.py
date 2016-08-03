from openFile import openFile

# This represents a genetic sequence from an HMM model.
class HmmSequence(object):
    def __init__(self, line):
        l = line.split();
        self.sequence = l[0];
        self.hmms = {}
        self.addHmm(l)
#        self.hmmFrom = int(l[4]);
#        self.hmmTo = int(l[5]);
        # Commenting these out to save space.
#        self.alignFrom = int(l[6]);
#        self.alignTo = int(l[7]);
        self.sequenceLength = int(l[10]);
        self.direction = l[11];

        self.alignedLength = abs(int(l[7]) - int(l[6])) # l[6/7] are alignFrom and alignTo. 
    
    def __str__(self):
        return self.hmmName;
    
    def isSplit(self, percent = .8):
        return float(self.alignedLength)/float(self.sequenceLength) <= percent

    # Tries to get the sequence's hmm.
    # Returns false if it doesn't exist.
# I don't use anymone.
#    def getHmm(self):
#        try:
#            return self.hammerTable.getHmm(self.hmmName)
#        except NameError: 
#            return False
    
    # Add the hmm to the list of hmms 
    # l is the split line from a hammer file.
    # l[2] = hmmName
    def addHmm(self, l):
        self.hmms[l[2]] = SequenceHmm(l)

    def getHmms(self):
        return set(self.hmms.keys())

    def getHmm(self, hmm):
        return self.hmms.get(hmm, None)
    # I don't think this is used. Will remove eventually.
#    @staticmethod
#    def createHmmSequences(fileLoc):
#        sequenceList = [];
#        with openFile(fileLoc) as sequences:
#            for line in sequences:
#                if line[0] != "#":
#                    sequenceList.append(HmmSequence(line));
#        return sequenceList;

class SequenceHmm(object):
    def __init__(self, l):
        self.hmmName = l[2]
        self.hmmFrom = l[4]
        self.hmmTo = l[5]
        self.alignFrom = l[6]
        self.alignTo = l[7]

