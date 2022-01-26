# Measures-independent-of-a-Model-of-Coding-DNA
My scripts to perform computations based on "Measures independent of a Model of Coding DNA" for my Masters project at Roderic Guigó's Lab in CRG.

1. Position Asymetry -> Assymetry.py
Measure how asymmetric is the distribution of nucleotides at the three triplet positions in a given window sequence.

2. Periodic Asymmetry Index -> Periodic_Assymetric_Index.py
Measure probability of finding pairs of the same nucleotide by comparing nucleotide i to nucleotide j (range i+1 to the end d of the sequence). Then computes homodimers per triplet position in each window.

3. Average Mutual Information -> Average_Mutual_Information.py
Measure the probability in the sequence of the pair of nucleotides i and j at a distance of k nucleotides.

4. Fourier Spectrum -> Fourier.c
Measure partial spectrum of a DNA sequence.

