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

sitei = 0 # site i
sitej = 3 # site j
scale = 4.0 # this number has to be such that max(|eigenvalues|)/scale < 1
npol = n # number of polynomials, energy resolution goes as 1/npol
ne = 300 # number of energies to calculate (between -scale and scale)
# returns energies and dos
(x1,y1r,y1i) = kpm.correlator0d(m,i=sitei,j=sitej,scale=scale,npol=npol,ne=ne) 


plt.plot(x1,y1r,label="KPM, real",c="red") # plot this dos
plt.plot(x1,y1i,label="KPM, imag",c="blue") # plot this dos

# use matrix inversion
iden = np.matrix(np.identity(n))
x2 = np.linspace(-3.,3.,100) # energies
delta = .02 # smearing
y2 = [((m.todense() - (x+1j*delta)*iden).I)[sitei,sitej] for x in x2]
y2 = np.array(y2) # in array form
y2r,y2i = -y2.real,y2.imag # real and imaginary parts
plt.scatter(x2,y2r,label="GF, real",c="red")
plt.scatter(x2,y2i,label="GF, imag",c="blue")

plt.title("Correlation function")
plt.legend()

plt.show()
