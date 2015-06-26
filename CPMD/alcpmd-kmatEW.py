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
# grep "out" from files
refOut = filter(lambda x:re.search(r'out$', x), refDir)
if len(refOut)>1:
  sys.exit("more then one output file in " + sys.argv[1])
# cat dir and out file for complete path
refPath = sys.argv[1] + "/" + refOut[0]
refFile = list(open(refPath,"r"))
# read file line by line into array
refLines = map(lambda each:each.strip("\n"), refFile)

# extract number of electrons
strNer = filter(lambda x:re.search(r'NUMBER\sOF\sELECTRONS', x), refLines)
Ner = [x.split(':')[1] for x in strNer][0]
Ner = int(float(re.findall("\d+.\d+", Ner)[0]))

# extract KS matrix dimension by flag 'KSDIMENSION'
refD = filter(lambda x:re.search(r'KSDIMENSION', x), refLines)
Dr = [int(s) for s in refD[0].split() if s.isdigit()][0]
#OCC = Ner/2
#if len(sys.argv) == 3:
#  VIR = Dr - OCC
#elif len(sys.argv) == 4:
#  Nvir = int(float(sys.argv[3]))
#  if Nvir <= Dr:
#    VIR = Nvir - OCC
#  else:
#    sys.exit("input dimension too big")
#else:
#  sys.exit("too many input arguments")

# extract matrix data by flag 'KSOUT'
refMat = filter(lambda x:re.search(r'KSOUT', x), refLines)
# remove 'KSOUT' from Mat varialbe
strKSr = [x.split('T ')[1] for x in refMat]


for i in range(0,len(strKSr)):
  x = strKSr[i]
  ksR = float(x.split(',')[0].translate(None, ''.join('(')))
  ksI = float(x.split(',')[1].translate(None, ''.join(')')))
  if (i==0):
    ksX = np.array([ksR + (ksI)*1j])
  else:
    ksX = np.vstack([ksX,ksR + (ksI)*1j])
    

#for x in strKSr:
#  ksR = float(x.split(',')[0].translate(None, ''.join('(')))
#  ksI = float(x.split(',')[1].translate(None, ''.join(')')))
##  if len(ksX==1):
#  ksX = ksR + (ksI)*1j
##  else:
#    
#  print ksX


#print strKSr[1].translate(None, ''.join(['(',')']))




# convert string array to np matrix of floats
#KSr = np.mat(np.array([float(x) for x in strKSr]).reshape((Dr,Dr)))
KSr = ksX.reshape(Dr,Dr)

Wr, Ur = np.linalg.eig(np.real(KSr))

Wr2 = np.diag(KSr)

print len(Wr), len(Wr2)
print KSr[0,0]*27.211396132
print Wr*27.211396132
print Wr2*27.211396132

#######################
##### data processing #
#######################
####
####
##### only closed shell system are considered
##### diagonalization of reference KS Hamiltonian
####Wr, Ur = np.linalg.eig(np.real(KSr))
####
##### rotate target KS Hamiltonian to reference MO
####Ht = Ur.T * KSt * Ur
####Wii = np.diag(Ht)
####Wia = Ht[0:OCC,OCC:VIR+OCC]
##### first order alchemical derivative
####d1E = 2*sum(Wii[0:OCC]-Wr[0:OCC]) + ESt - ESr
##### second order alchemical derivative
####d2E = 0
####for i in range(0,OCC):
####  for a in range(0,VIR):
####    #d2E=0
####    d2E += 2*Wia[i,a]**2/(Wr[i] - Wr[OCC+a])
####
####
#####print "KSr:"
#####print KSr
#####print "Ur"
#####print Ur
#####print "Wr"
#####print Wr
####
#####print Et, d1E, d2E
####print Et+d1E, Et+d1E+d2E
#####print Et+d1E+d2E
#####print d1E.real, d2E.real
#####print d1E
#####print d2E
