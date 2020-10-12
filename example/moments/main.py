import sys
sys.path.append("../../src") # path to KPM library
import kpm  # library with the kernel polynomial method
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix # class for sparse matrices
import numpy as np
import time


###############################################
## this a test with a 1d tight binding model ##
###############################################
n = 1000 # dimension of the matrix
print("Dimension of the matrix",n)
ii = np.arange(0,n-2) # indexes from 0 to n-2
cols = np.concatenate([ii,ii+1]) # index for rows
rows = np.concatenate([ii+1,ii]) # index for cols
data = np.zeros(cols.shape[0]) + 1. # tight binding parameter

m = csc_matrix((data,(rows,cols)),shape=(n,n)) # create the sparse matrix

site = n//2 # the entry where you want the dos, n/2 is an atom in the middle
scale = 4.0 # this number has to be such that max(|eigenvalues|)/scale < 1
npol = 300 # number of polynomials, energy resolution goes as 1/npol
ne = 300 # number of energies to calculate (between -scale and scale)
# returns energies and dos

# this function computes the different moments using the KPM
mus = kpm.get_moments_ldos(m,i=site,n=npol,scale=scale)

# once you have the moments, you can reconstruct the DOS
# using the Chebyshev relations
x = np.linspace(-0.9,0.9,400)
y = kpm.generate_profile(mus.real,x,kernel="jackson").real


# plot the result
plt.plot(x*scale,y/scale) # plot this dos
plt.show()
