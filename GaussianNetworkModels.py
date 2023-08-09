# Gaussian Network Model (GNM)
# It gives the magnitude of the motions, which means it does not contain the direction information.
# Most of the contribution of motion is stemmed from the lowest eigenvalue.

from prody import *
from matplotlib.pylab import *
import py3Dmol
import pandas as pd
import numpy as np
# Enter the PDB name and selection as x= PDBID and y= protein.
protein = parsePDB("input.pdb").select("protein")
# Select Calpha atoms within the PDB.
protein_calpha = protein.select("calpha")
# Build your GNM.
gnm = GNM('protein_calpha')
# Generate your Kirchhof matirx.
gnm.buildKirchhoff(protein_calpha, gamma=1, cutoff=7.3)
# Instead of n, enter the number of mode.
gnm.calcModes(len(protein_calpha))
print(gnm.getEigvals())
eigVals = gnm.getEigvals()
eigVecs = gnm.getEigvecs()
collectivity = []
i = 0
collectivity_ls = list()
while i <= len(eigVecs) - 1 :
    collectivity = (1/len(protein_calpha)) * np.exp(-np.sum(np.square(eigVecs[i,:])*np.log(np.square(eigVecs[i,:]))))
    collectivity_ls.append(collectivity)
    i = i + 1
#print(collectivity_ls)
savetxt('eigVals.txt', eigVals,  delimiter=',')
savetxt('eigVecs.txt', eigVecs, delimiter = ',')
eVec = open('eigVecs.txt', 'r')
resnums = protein_calpha.getResnums()
print('The contribution rates of the eigenvalues...')
contributionOfeigVals = 1/eigVals
print(contributionOfeigVals)
i = 0
print('Uncovering the hinge residues...')
while i <= 2 :
  hinges=calcHinges(gnm[i])
  print(resnums[hinges]+127)
  i = i + 1
  
fig,ax=plt.subplots()
fig.canvas.draw()
labels = [item.get_text() for item in ax.get_xticklabels()]

newLabels = [127, 177, 227, 277, 327, 377, 427]
i=0

while i <= len(newLabels)-1 :

    if i <= len(labels)-1:
        print(i)
        labels[i] = newLabels[i]

    else:
        labels.append(newLabels[i])

    i=i+1
ax.set_xticklabels(labels)
j = 0
system = ['Mode 1', 'Mode 2', 'Mode 3']
color = ['b', 'r', 'g']
while j <= 2 :
    showSqFlucts(gnm[j], hinges=True, label=system[j])
    j = j + 1
plt.title("")
plt.xlabel("Residue index", fontweight="bold")
plt.ylabel("Square fluctuations (Å$^2$)", fontweight="bold")
plt.legend(loc=1)
plt.show()
plt.savefig('GNMsqrFlucFirstCluster.png', dpi=300)
