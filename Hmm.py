import sys
import numpy
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as pyplot
from openFile import openFile

# Represents an HMM model. 
# Can come from 2 sources
#    1) A skew file, which hold all the data
#    2) Just the name, and then build up data by adding hammer table lines for data.
class Hmm(object):
    def __init__(self, line=None):
        l = line.strip().split(",");
        if len(l) == 11: # If the model is being built from a skew file.
            self.model = l[0];
            self.modelLength = int(l[1]);
            self.snumpyCoverageCoefficent = l[2];
            self.hits = l[3];
            self.percentCoverage = l[4];
            self.avgerageCoverageDepth = l[5];
            self.shannonEntropy = l[6];
            self.metricEntropy = l[7];
            self.minMax = l[8];
            self.l2normDeviation = l[9];
            self.vector = [int(x) for x in l[10].split(" ")];
        elif len(l) == 1: # If only the model name is passed in. 
            self.model = l[0];
            self.last = 0;
            self.vector = []

    # Takes a line from a hammer table file and uses it to build coverage 
    # data for the model.
    def addLine(self, l):
        if l[2] != self.model:
            raise Exception("Line does not go to this model. Expected: " + self.model + " but got: " + l[2])
        
        if int(l[5]) > self.last:
            self.vector.extend([0]*(int(l[5]) - self.last))
            self.last = int(l[5])
            self.modelLength = len(self.vector)
        
        for x in range(int(l[4]), int(l[5]) + 1):
            self.vector[x-1] += 1

    # Calculates the coverage, and figures out where areas of low covarge are.
    def calcBlocks(self, cutoff=0.75):
        self.breakpoints = [0] * self.modelLength
        self.coverageMean = numpy.average(self.vector)
        self.coverageStd = numpy.std(self.vector)
        self.missingBlocks = [0] * len(self.vector)
        for i, value in enumerate(self.vector):
            if value < self.coverageMean * cutoff:
                self.missingBlocks[i] = 100
            else:
                self.missingBlocks[i] = 0

    # Breakpoints are where sequences start and stop. We thought this would be
    # helpful, but it was not. Probably won't be used.
    def addBreakpoints(self, start, end):
        self.breakpoints[start-1] += 1;
        self.breakpoints[end-1] += 1;

    # Creates and saves a graph that shows the coverage and areas of low 
    # coverage. 
    def drawVectorAndBreakpoints(self, location = "graphs/"):
        # Testing to see if this works. It may not.
        if self.coverageMean == 0:
            return
        pyplot.plot(self.vector);
        pyplot.plot(self.breakpoints);
        pyplot.plot([x * self.coverageMean * .8 for x in self.missingBlocks]);
        pyplot.axhline(self.coverageMean);
        pyplot.axhline(self.coverageMean*1.1);
        pyplot.axhline(self.coverageMean*.75);
        pyplot.xlabel("Base");
        pyplot.ylabel("Count");
        pyplot.title(self.model + ": " + str(self.modelLength) 
                + " mean: " + str(self.coverageMean)
                + " std: " + str(self.coverageStd))
        pyplot.draw();
        pyplot.savefig(location + "coverage-and-breakpoints-" + self.model + ".png");
        pyplot.close()

    def __str__(self):
        return self.model;
    def __repr__(self):
        return self.__str__(); 

    # Create a dictionary with HMMs from a skew file.
    @staticmethod
    def createHmms(fileLoc):
        hmmDic = {};
        with openFile(fileLoc) as hmms:
            next(hmms); # Skip the header
            for line in hmms:
                newHmm = Hmm(line);
                newHmm.calcBlocks();
                hmmDic[newHmm.model] = newHmm;
        return hmmDic;

    # Creates and dictionary with HMMs from a hammer file.
    # It also does the calculations needed to create graphs.
    @staticmethod
    def createHmmsFromHammer(fileLoc):
        hmmDic = {}
        with openFile(fileLoc) as hammer:
            # Skip the first two lines.
            next(hammer)
            next(hammer)
            for line in hammer:
                if line.strip()[0] == "#": # Skip comment lines.
                    continue
                l = line.strip().split();
                hmmName = l[2]

                if not hmmDic.has_key(hmmName):
                    hmmDic[hmmName] = Hmm(hmmName)

                hmmDic[hmmName].addLine(l)
        for hmm in hmmDic:
            hmmDic[hmm].calcBlocks()
        return hmmDic;

    # I don't think I used this any more. Will remove eventually.
    @staticmethod
    def createHmmSequences(fileLoc):
        hmm = {}
        with openFile(fileLoc) as hmmFile:
            for line in hmmFile:
                l = line.split(" ");
                hmm[l[0]] = l[1:]
        return hmm



def main():
   Hmm.createHmmsFromHammer(sys.argv[1])

if __name__ == "__main__":
    main()
