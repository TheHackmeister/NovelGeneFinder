import sys
import numpy
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as pyplot
from Hmm import Hmm
from HmmSequence import HmmSequence

def calcDirection(fileLoc):
    plus = 0;
    minus = 0;
    other = 0;
    wrong = 0;
    with open(fileLoc) as hammer:
        next(hammer);
        next(hammer);

        for line in hammer:
            try: 
                if line.split()[11] == "-":
                    minus += 1;
                elif line.split()[11] == "+":
                    plus += 1;
                else:
                    other += 1;
            except IndexError:
                wrong +=1;

    return plus, minus, other, wrong;

def createHmmSequenceAlignementCount(hmmSequences):
    pyplot.hist([x.alignedLength for x in hmmSequences], bins=numpy.arange(0,102));
    pyplot.xlabel("Number of bases aligned to HMM model");
    pyplot.ylabel("Count");
    pyplot.draw();
    pyplot.savefig("alignmentCount.png");


def pickOeaSplitReads(reads, percent = .8):
    oeaReads = []
    for sequence in reads:
        if sequence.isSplit(percent):
            oeaReads.append(sequence)

    return oeaReads

def main():
    hmms = Hmm.createHmmsFromHammer(sys.argv[2]); 
    hmmSequences = HmmSequence.createHmmSequences(sys.argv[2]);


    notFound = {}
    for sequence in hmmSequences:
        if sequence.isSplit():
            try:
                hmms[sequence.hmmName].addBreakpoints(sequence.hmmFrom, sequence.hmmTo)
            except KeyError:
                notFound[sequence.hmmName] = True

    print notFound

    for hmm in hmms:
        hmms[hmm].drawVectorAndBreakpoints();


if __name__ == "__main__":
    main()

# On the 9.tblout.scan I get the following.
# (462207, 53540, 1, 9)
# That is, + is responsible for 89% of the reads.
# print calcDirection(sys.argv[1]);



