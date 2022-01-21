# computes entropy of DNA sequences
# see Fickett and Tung, NAR 20:6441-6450 (1992)
# USAGE=python3 Assymetry.py [wlen] [wstep] [fasta]
import sys
from typing import Sequence

def assymetry(wseq, wlen):
    # get base composition at codon positions
    nuclwind2count = {}
    i = 1
    for nucl in wseq:
        if (nucl, i%3) not in nuclwind2count.keys():
            nuclwind2count[(nucl, i%3)] = 1
        else:
            nuclwind2count[(nucl, i%3)] += 1
        i += 1
    # get probabilities. dividing by l/3
    for key in nuclwind2count:
        nuclwind2count[key]/= wlen/3
    # compute assymetries
    # average nucleotide frequency per position
    anf = {}
    for nucl in nucl2prob:
        for i in range(0,3):
            if nucl not in anf.keys():
                anf[nucl] = 0
            try:
                anf[nucl] += nuclwind2count[(nucl, i)]
            except:
                pass
        anf[nucl] /= 3
    # now assymetries
    assym = {}
    for nucl in nucl2prob:
        for i in range(0, 3):
            if nucl not in assym.keys():
                assym[nucl] = 0.0
            if (nucl, i) in nuclwind2count.keys():
                assym[nucl] += (nuclwind2count[(nucl, i)]-anf[nucl])**2
            
    return (assym["A"]+assym["C"]+assym["G"]+assym["T"])


global nucl2prob
nucl2prob = {'A': 1, 'T':1, 'C':1, 'G':1}
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
    print(i+1, i+WLEN, round(assymetry(wseq, WLEN), 7))
    i += WSTEP

# HUMHBB.2ex 
# GCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGG
