import sys

def varifyMatchingPairs(fileLoc):
    with open(fileLoc) as reads:
        for line in reads:
            read1 = line.split(" ");
            next(reads);
            read2 = next(reads).split(" ");
            next(reads);

            if read1[0] != read2[0]:
                raise Exception("Missmatch!");

            if read1[1].split(":")[0] != "1":
                raise Exception("Read 1 problem!");

            if read2[1].split(":")[0] != "2":
                raise Exception("Read 2 problem!");

def verifySequenceLength(fileLoc):
    with open(fileLoc) as reads:
        for line in reads:
            line = next(reads);
            try:
                if len(line) != 102:
                    raise Exception("No go!" + str(len(line)));
            except NameError:
                    length = len(line);


varifyMatchingPairs(sys.argv[1]);
print "Finished matching pairs";
verifySequenceLength(sys.argv[1]);
print "Finished sequence length";



print "Done!"


