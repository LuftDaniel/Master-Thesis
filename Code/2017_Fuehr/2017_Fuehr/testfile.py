import sovi_bib as sovi
import numpy as np


grad = [1.0]
defo = [3.0, 2.0, 1.0]

zero = np.zeros([3,3])

mem = sovi.bfgs_memory(zero, zero, 3)

for i in range(2):
    i = i+1
    print(i)