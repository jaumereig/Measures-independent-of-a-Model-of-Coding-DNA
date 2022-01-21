import sys
import numpy as np
from numpy.fft import *

def fourier_genome(wseq, wlen):
    # Array for each nucleotide of 0 with size wseq-1
    xA=[0 for i in range(0,len(wseq))]
    xT=[0 for i in range(0,len(wseq))]
    xG=[0 for i in range(0,len(wseq))]
    xC=[0 for i in range(0,len(wseq))]

    # If position i in each array corresponds to its nucleotide write 1
    for i in range(0, wlen):
        if wseq[i]=="A":xA[i]=1
        if wseq[i]=="T":xT[i]=1
        if wseq[i]=="G":xG[i]=1
        if wseq[i]=="C":xC[i]=1
     
    # This function computes the one-dimensional *n*-point discrete Fourier Transform (DFT) with the efficient 
    # Fast Fourier Transform (FFT) algorithm [CT]. 
    xhatA=abs(fft(xA)) # return absolute value of argument
    xhatG=abs(fft(xG)) # output is in array of absolute numbers corresponding to fft() in each position
    xhatT=abs(fft(xT))
    xhatC=abs(fft(xC))

    l = int(len(xhatA)) # i.e. 223 when using HUMHBB.2ex.fa
    # Now we do sum of squares of each position per each nucleotide in order to have a single array resulting as a combination  of each nucl
    xhat=[(xhatA[i]**2+xhatT[i]**2+xhatC[i]**2+xhatG[i]**2)*2/len(wseq) for i in range(int(l/2),l)]
    #return(xhat,[float(i)/len(xhat) for i in range(0,len(xhat))])
    xhat = [float(i)/len(xhat) for i in range(0,len(xhat))]
    return xhat

# def P(spectrum):
#     x=len(spectrum)
#     k=int(min(x/3,50))
#     peak=max(spectrum[int(x/3)-k:int(x/3)+k])
#     p=x*peak/sum(spectrum)
#     return(p)


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
    #ff_p, freq_p = fourier_genome(wseq, WLEN)
    print(i+1, i+WLEN, fourier_genome(wseq, WLEN))
    i += WSTEP