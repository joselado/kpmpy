import sys
sys.path.append("../src") # path to KPM library
import kpm  # library with the kernel polynomial method
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix # class for sparse matrices
import numpy as np


###############################################
## this a test with a 1d tight binding model ##
###############################################
compare = False # compare with the normal method
n = 200 # dimension of the matrix
ii = np.arange(0,n-2) # indexes from 0 to n-2
cols = np.concatenate([ii,ii+1]) # index for rows
rows = np.concatenate([ii+1,ii]) # index for cols
data = np.zeros(cols.shape[0]) + 1. # tight binding parameter

m = csc_matrix((data,(rows,cols)),shape=(n,n)) # create the sparse matrix

site = n//2 # the entry where you want the dos, n/2 is an atom in the middle
scale = 4.0 # this number has to be such that max(|eigenvalues|)/scale < 1
npol = 100 # number of polynomials, energy resolution goes as 1/npol
ne = 300 # number of energies to calculate (between -scale and scale)
# returns energies and dos
(x1,y1) = kpm.ldos0d(m,i=site,scale=scale,npol=npol,ne=ne) 


plt.plot(x1,y1,label="KPM") # plot this dos

# use matrix inversion
iden = np.matrix(np.identity(n))
x2 = np.linspace(-3.,3.,400) # energies
delta = .05 # smearing
y2 = [((m.todense() - (x+1j*delta)*iden).I)[site,site].imag for x in x2]
y2 = np.array(y2)/np.pi # pi factor
plt.scatter(x2,y2,label="Green function")



plt.show()
