import sys
import os
import linecache
from itertools import izip, izip_longest



# Build hmm as a dictionary with the sequence as the key. 
# This provides O(1) lookup. 
hmm = {} 
with open(sys.argv[1]) as hmmFile: 
    for line in hmmFile:
        split = line.split(" ");
        hmm[split[0]] = split[1:];

flag = False;

fOrphanFile = sys.argv[2] + ".orphans";
rOrphanFile = sys.argv[3] + ".orphans";
if os.path.exists(fOrphanFile):
    raise Exception("The forward orphan file already exists.");
if os.path.exists(rOrphanFile):
    raise Exception("The reverse orphan file already exists.");

with open(sys.argv[2]) as fastqf, open(sys.argv[3]) as fastqr, open(fOrphanFile, "w+") as fOrphan, open(rOrphanFile, "w+") as rOrphan:
    for fLine, rLine in izip_longest(fastqf, fastqr):
        if(fLine.startswith("@")):
            if(fLine.split(" ")[0] != rLine.split(" ")[0]):
                raise Exception("The fastq file contains reads in mismatching orders.");
            fSequence = next(fastqf); 
            rSequence = next(fastqr); 
            if not hmm.get(fSequence[:-1], False) and not hmm.get(rSequence[:-1], False): # last char is a new line. 
                fOrphan.write(fLine + fSequence + next(fastqf) + next(fastqf));
                rOrphan.write(rLine + rSequence + next(fastqr) + next(fastqr));

print "Done";
