# computes Periodic Assymetry Index (Konopka, 1990)
# USAGE= python3 -f Periodic_Assymetric_Index.py window_length window_step sequence_file
import sys

def PAI(wseq, wlen):
    # relation between 1 bp i and consecutive bps j
    pairnucldis2count = {}
    # dict stores nucleotide pairs while i fixed and distances between them, then counts occurrences of each
    for i in range(1, wlen):
        j = i+1
        n1 = wseq[i-1]
        while j <= wlen:
            if (n1, wseq[j-1], j-i-1) not in pairnucldis2count.keys():
                pairnucldis2count[(n1, wseq[j-1], j-i-1)] = 1
            else:
                pairnucldis2count[(n1, wseq[j-1], j-i-1)] += 1
            j += 1

    # compute number of homodimers at k%3 = 0, 1, 2 (codon positions) per window
    A3 = {}
    for start in range(3):
        A3[start] = 0
    for i in nucleotides:
        for k in range(wlen-1):
            try:
                A3[k%3] += pairnucldis2count[(i, i, k)]
            except:
                pass
    
    MA3 = max(max(A3[0], A3[1]), A3[2]) / (min(min(A3[0], A3[1]), A3[2]) +1)

    return MA3
        
global nucleotides
nucleotides = ['A', 'C', 'G', 'T']
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
    print(i+1, i+WLEN, round(PAI(wseq, WLEN), 7))
    i += WSTEP