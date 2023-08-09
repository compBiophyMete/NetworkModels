# Anisotropic network model

from prody import *
from matplotlib.pylab import *
import py3Dmol
import pandas as pd
import numpy as np
# Enter the PDB name and selection as x= PDBID and y= protein.
print('Parsing the PDB file...')
protein = parsePDB("cluster1_allAtom.pdb").select("protein")
# Select Calpha atoms within the PDB.
print('Selecting the Calpha atoms...')
protein_calpha = protein.select("calpha")
# Build your ANM.
print('Building Heissian matrix...')
anm = ANM('protein_calpha')
anm.buildHessian(protein_calpha)

# Instead of n, enter the number of mode.
print('Calculating the eigenvalues and eigenvectors...')
anm.calcModes(len(protein_calpha))
print(anm.getEigvals())
eigVals = anm.getEigvals()
print('The number of eigenvalues:',len(eigVals))
eigVecs = anm.getEigvecs()
collectivity = []
i = 0
#Perturbation Response Scanning
print('Applying perturbations...')
prs_mat, eff, sens= calcPerturbResponse(anm[0:2])
#30, 67, 231,61
TriadIndeces = [67, 231]
for TriadIndex in TriadIndeces :
    np.savetxt('out' + str(TriadIndex+127) + '.txt',prs_mat[:,TriadIndex-1])
