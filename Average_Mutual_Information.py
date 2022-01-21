# computes Average Mutual Information. Simplified way. Ivo Grosse
# USAGE= gawk -f AMI.awk window_length window_step sequence_file
import sys
import numpy as np

def AMI(wseq, wlen):
    # get base composition at codon positions
    nucl2numpos = {} # number of each nucleotide per position in triplet
    nucl2numtot = {} # number of each nucleotide in total
    for n in nucleotides:
        for i in range(3):
            nucl2numpos[(n, i)] = 0
        nucl2numtot[n] = 0
    i = 0
    for nucl in wseq:
        # if (nucl, i%3) not in nucl2numpos.keys():
        #     nucl2numpos[(nucl, i%3)] = 1
        # else:
        nucl2numpos[(nucl, i%3)] += 1
        i += 1
        nucl2numtot[nucl] += 1
    # for nucl in wseq:
    #     # if (nucl) not in nucl2numtot.keys():
    #     #     nucl2numtot[nucl] = 1
    #     # else:
    #     nucl2numtot[nucl] += 1
    #     #i += 1

    # # get probabilities
    for key in nucl2numpos:
        nucl2numpos[key]/= wlen/3 # frequencies per codon position of each nucleotide
    for key in nucl2numtot:
        nucl2numtot[key]/=wlen # frequencies of nucleotide along all window

    # Now compute dinucleotide probabilities  (Herzel and Grosse)
    # inframe
    dinuc2prob_in = {}
    for i in nucleotides:
        for j in nucleotides:
            # for k in range(3):
            #     if (nucleotides[i],k) not in nucl2numpos.keys():
            #         nucl2numpos[(nucleotides[i],k)] = 0

            #     if (nucleotides[j],k) not in nucl2numpos.keys():
            #         nucl2numpos[(nucleotides[j],k)] = 0

            dinuc2prob_in[(i, j)] = ((nucl2numpos[(i,0)]*nucl2numpos[(j,0)])+(nucl2numpos[(i,1)]*nucl2numpos[(j,1)])+(nucl2numpos[(i,2)]*nucl2numpos[(j,2)])) / 3




    # outframe
    dinuc2prob_out = {}
    for i in nucleotides:
        for j in nucleotides:
            dinuc2prob_out[(i,j)] = ((nucl2numpos[(i,0)]*nucl2numpos[(j,1)])+(nucl2numpos[(i,1)]*nucl2numpos[(j,2)])+(nucl2numpos[(i,2)]*nucl2numpos[(j,0)])) / 3
                # F_out[b[i],b[j]] = ((G[b[i],0]*G[b[j],2])+(G[b[i],1]*G[b[j],0])+(G[b[i],2]*G[b[j],1])) / 3
                
    # Compute I_in and I_out (Average Mutual Information)
    ami2prob_in = 0
    for i in nucleotides:
        for j in nucleotides:
            ami2prob_in += (dinuc2prob_in[(i,j)]*np.log2(dinuc2prob_in[(i,j)]/(nucl2numtot[i]*nucl2numtot[j])))
    ami2prob_out = 0
    for i in nucleotides:
        for j in nucleotides:
            ami2prob_out += (dinuc2prob_out[(i,j)]*np.log2(dinuc2prob_out[(i,j)]/(nucl2numtot[i]*nucl2numtot[j])))

    return (2/3)*ami2prob_out + (1/3) * ami2prob_in
    #return (ami2prob_in/ami2prob_out)


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
    print(i+1, i+WLEN, AMI(wseq, WLEN))
    i += WSTEP