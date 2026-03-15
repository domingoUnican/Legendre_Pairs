from itertools import product
import numpy as np

alphabet = [-1,1] # This are the elements which are made the sequences
# Example of shifts of a sequence (-1,-1,-1,1)
# sequence shift by one -> (-1, -1, 1, -1)
# sequence shift by two -> (-1, 1, -1, -1)
# sequence shift by three -> (1, -1, -1, -1)
def PAF_i(seq, i):
    """
    The PAF function takes a sequence of length n and returns a sequence of length n.
    This function returns the element i of the PAF
    """
    return vector(seq).dot_product(vector(seq[i:]+ seq[:i]))
#examples of calculation of PAF((-1,-1,-1,1),i)
# if i==0 -> 4
# if i==1 -> 0
# if i==2 -> 0
# if i==3 -> 0

def PAF(seq):
    """
    This is fast version using the Fast Fourier Transform, 
    but the first element is remove because it is always len(seq).
    """
    temp = np.fft.fft(seq)
    temp *= np.roll(temp[::-1],1) # This is the formula for the PAF
    temp = np.fft.ifft(temp)
    return tuple(np.array(np.rint(temp),int)[1:])
# This version uses the FFT and it has complexity O(n log n)
#which compares with PAF_i, which has complexity O(n) for each call,
# O(n^2).

def Compresion(seq, m):
    """
    This is function that compress a sequence, summing the elements
    """
    L = len(seq)
    result = []
    for i in range(m):
        result.append(sum((seq[i+m*j] for j in range(L//m))))
    return tuple(result)

# Compresion((-1,-1,-1,1),2)
# (-2,0, -2,0)
# Because it is periodic (-2,0)
# Compresion((-1,-1,-1,1,1,1),3)
# (0,0,0)
# Compresion((-1,-1,-1,1,1,1),2)
# (-1,1)
dft_set = dict()
pairs = set()
p = 5
L = 9
for i in product(range(5), repeat=L):
    dft_set[PAF(i)] = i
    temp = tuple([p*(p*L-3)//2-j for j in PAF(i)])
    if temp in dft_set:
        pairs.add((i, dft_set[temp]))
g = open("pairs_9.txt", "w")
for pair in pairs:
    a, b = pair
    g.write(f"{' '.join([str(i) for i in a])}\n{' '.join([str(i) for i in b])}\n")
g.close()
