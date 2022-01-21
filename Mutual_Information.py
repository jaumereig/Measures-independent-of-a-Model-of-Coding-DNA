
import sys
import numpy as np




global nucleotides
nucleotides = ['A', 'C', 'G', 'T']
KMAX=50
WLEN = int(sys.argv[1])
WSTEP = int(sys.argv[2])
sequence_file = sys.argv[3]

sequence =''
with open(sequence_file, 'r') as sf:
    for line in sf.readlines():
        if line[0] != '>':
            sequence += line[:-1]
    # for line in sf:
    #     sequence += line.split()[1]


lseq = len(sequence)
lwseq=lseq-WLEN
i = 0
while i <= lwseq:
    wseq = sequence[i:i+WLEN]
    print(i+1, i+WLEN, AMI(wseq, WLEN))
    i += WSTEP