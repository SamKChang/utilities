#!/usr/bin/python
#
# comput Coulomb matrix from xyz file
# input: foo.xyz baseDimension
# output: reshaped column vector of upper triangular Coulomb matrix
#
# output vector can be used to evaluate norm

def coulombM(fXYZ,baseDim):
  import numpy as np
  XYZ = open(fXYZ,"r")
  na = int(XYZ.readline())
  XYZ.readline()
  
  def atom_a2n(ATOM):
    return {
      'H':1,
      'C':6,
      'N':7,
      'O':8,
      'Al':13,
      'S':16,
      'Ga':31,
      'As':33,
    } [ATOM]

  # read xyz file to data
  for i in range(0,na):
    line=XYZ.readline().split()
    if(i==0):
      data=np.array(
        [atom_a2n(line[0]),line[1],line[2],line[3]]
      ).astype(float)
    else:
      tmp=np.array(
        [atom_a2n(line[0]),line[1],line[2],line[3]]
      ).astype(float)
      data=np.vstack([data,tmp])
  
  # construct upper triangular Coulomb matrix
  CM=np.zeros((na,na))
  for i in range(0,na):
    for j in range(i+1,na):
      Rij=np.linalg.norm(data[i][1:]-data[j][1:])
      CM[i][j]=data[i][0]*data[j][0]/Rij
  # full Coulomb matrix
  CM = CM + np.transpose(CM)
  for i in range(0,na):
    CM[i][i]=0.5*data[i][0]**2.4
  
  # row norm of Coulomb matrix
  NORM=np.zeros((na,2))
  for i in range(0,na):
    NORM[i,0]=np.linalg.norm(CM[i])
    NORM[i,1]=i
  
  # sort NORM by row norm
  NORM=NORM[NORM[:,0].argsort()]
  # reverse array order
  NORM=NORM[::-1]
  # extract new row order
  sortRow=NORM[:,1].astype(int)
  # rearrange Coulomb matrix
  CM=CM[:,sortRow][sortRow]
  
  # fill Coulomb matrix with zeros
  dim=int(baseDim)
  d=dim-na
  B=np.zeros((na,d))
  C=np.zeros((d,d))
  D=np.hstack([np.transpose(B),C])
  E=np.hstack([CM,B])
  CM=np.vstack([E,D])
  
  #print CM
  dimV=(dim**2 - dim)/2 + dim
  CV=np.zeros(dimV)
  s=0
  for i in range(0,dim):
    for j in range(i,dim):
      CV[s]=CM[i][j]
      s=s+1

  return CV

if __name__ == "__main__":
  import sys
  CV=coulombM(sys.argv[1],sys.argv[2])
  dim=int(sys.argv[2])
  dimV=(dim**2 - dim)/2 + dim
  for i in range(0,dimV):
    print CV[i]
