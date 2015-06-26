#!/usr/bin/python

# calculate 1st and 2nd order corrections from
# output KS Hamiltonian with FLAG KSOUT, KSDIMENSION

import sys                       # command line variable
import re                        # regular expression
import numpy as np               # numberical routines
from os import listdir           # handling files
from os.path import isfile, join # handling files

############################
# reading CPMD output file #
############################
# list files from directories specified by command line variables
refDir = [f for f in listdir(sys.argv[1]) if isfile(join(sys.argv[1],f))]
tarDir = [f for f in listdir(sys.argv[2]) if isfile(join(sys.argv[2],f))]
# grep "out" from files
refOut = filter(lambda x:re.search(r'out$', x), refDir)
tarOut = filter(lambda x:re.search(r'out$', x), tarDir)
if len(refOut)>1:
  sys.exit("more then one output file in " + sys.argv[1])
if len(tarOut)>1:
  sys.exit("more then one output file in " + sys.argv[2])
# cat dir and out file for complete path
refPath = sys.argv[1] + "/" + refOut[0]
tarPath = sys.argv[2] + "/" + tarOut[0]
refFile = list(open(refPath,"r"))
tarFile = list(open(tarPath,"r"))
# read file line by line into array
refLines = map(lambda each:each.strip("\n"), refFile)
tarLines = map(lambda each:each.strip("\n"), tarFile)

# extract number of electrons
strNer = filter(lambda x:re.search(r'NUMBER\sOF\sELECTRONS', x), refLines)
strNet = filter(lambda x:re.search(r'NUMBER\sOF\sELECTRONS', x), tarLines)
Ner = [x.split(':')[1] for x in strNer][0]
Net = [x.split(':')[1] for x in strNet][0]
Ner = int(float(re.findall("\d+.\d+", Ner)[0]))
Net = int(float(re.findall("\d+.\d+", Net)[0]))
if(Ner != Net) or (Ner%2 == 1):
  print "ref Ne:", Ner
  print "tar Ne:", Net
  sys.exit("only implemented for closed shell/isoelectron paths")

# extract electrostatic energy
#refES = filter(lambda x:re.search(r'ELECTROSTATIC', x), refLines)
#tarES = filter(lambda x:re.search(r'ELECTROSTATIC', x), tarLines)
#strESr = [x.split('Y = ')[1] for x in refES]
#strESt = [x.split('Y = ')[1] for x in tarES]
#ESr = float([x.split(' A.U.')[0] for x in strESr][1])
#ESt = float([x.split(' A.U.')[0] for x in strESt][1])

# extract reference total energy
#refEt = filter(lambda x:re.search(r'TOTAL ENERGY', x), refLines)
#strEt = [x.split('Y = ')[1] for x in refEt]
#Et = float([x.split(' A.U.')[0] for x in strEt][1])

# extract KS matrix dimension by flag 'KSDIMENSION'
refD = filter(lambda x:re.search(r'KSDIMENSION', x), refLines)
tarD = filter(lambda x:re.search(r'KSDIMENSION', x), tarLines)
Dr = [int(s) for s in refD[0].split() if s.isdigit()][0]
Dt = [int(s) for s in refD[0].split() if s.isdigit()][0]
if(Dr != Dt):
  print "refD:", Dr
  print "tarD:", Dt
  sys.exit("Kohn-Sham matrices have different dimention")
OCC = Ner/2
if len(sys.argv) == 3:
  VIR = Dr - OCC
elif len(sys.argv) == 4:
  Nvir = int(float(sys.argv[3]))
  if Nvir <= Dr:
    VIR = Nvir - OCC
  else:
    sys.exit("input dimension too big")
else:
  sys.exit("too many input arguments")

# extract matrix data by flag 'KSOUT'
refMat = filter(lambda x:re.search(r'KSOUT', x), refLines)
tarMat = filter(lambda x:re.search(r'KSOUT', x), tarLines)
# remove 'KSOUT' from Mat varialbe
strKSr = [x.split('T ')[1] for x in refMat]
strKSt = [x.split('T ')[1] for x in tarMat]
# convert string array to np matrix of floats
KSr = np.mat(np.array([float(x) for x in strKSr]).reshape((Dr,Dr)))
KSt = np.mat(np.array([float(x) for x in strKSt]).reshape((Dr,Dr)))

###################
# data processing #
###################
# only closed shell system are considered
# diagonalization of reference KS Hamiltonian
Wr, Ur = np.linalg.eig(KSr)
# rotate target KS Hamiltonian to reference MO
Ht = Ur.T * KSt * Ur
Wii = np.diag(Ht)
Wia = Ht[0:OCC,OCC:VIR+OCC]
# first order alchemical derivative
d1E = 2*sum(Wii[0:OCC]-Wr[0:OCC])# + ESt - ESr
# second order alchemical derivative
d2E = 0
for i in range(0,OCC):
  for a in range(0,VIR):
    d2E += 2*Wia[i,a]**2/float(Wr[i] - Wr[OCC+a])

#print Et, d1E, d2E
#print Et+d1E, Et+d1E+d2E
print d1E
