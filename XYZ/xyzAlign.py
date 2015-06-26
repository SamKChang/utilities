#!/usr/bin/python
# http://nghiaho.com/?page_id=671

from numpy import *
from math import sqrt
import re
import sys

aNameA = []
aCrdA = []
listA = []
Ca = array([0.0,0.0,0.0])
aNameB = []
aCrdB = []
listB = []
Cb = array([0.0,0.0,0.0])

###############
# print usage #
###############
def usage():
  print ""
  print "xyzAlign.py align A.xyz to B.xyz using least square approach." 
  print "At least 3 atoms from A.xyz and B.xyz are needed."
  print ""
  print "usage: xyzAlign.py A.xyz a1 a2 a3 ... B.xyz b1 b2 b3 ..."
  print ""
  sys.exit()

def printXYZ(NA,CRD,NAME):
  print NA
  print ""
  for i in range(0,NA):
    print("%-2s % 12.8f % 12.8f % 12.8f") % (NAME[i], CRD[i][0], CRD[i][1], CRD[i][2])
  
# check number of input
if len(sys.argv) < 9:
  usage()

# reading A.xyz 
xyzIn = open(sys.argv[1],'r')
Na = int(xyzIn.readline())
xyzIn.readline()
for i in range(0,Na):
  data = re.sub("[\n\t]", "",xyzIn.readline()).split(' ')
  # remove empty elements
  data = filter(None, data)
  aNameA.append(data[0])
  crd = [float(data[1]),float(data[2]),float(data[3])]
  aCrdA.append(crd)
#  Ca += array(crd)
aCrdA = vstack(aCrdA)
#Ca /= float(Na)

# reading aligning atoms of A.xyz
itr = 2
xyz = 0
while xyz==0:
  if sys.argv[itr].find('xyz') != -1:
    xyz += 1
  else:
    listA.append(int(sys.argv[itr])-1)
    itr += 1
if(len(listA)<3):
  usage()

# reading B.xyz
xyzIn = open(sys.argv[itr],'r')
Nb = int(xyzIn.readline())
xyzIn.readline()
for i in range(0,Nb):
  data = re.sub("[\n\t]", "",xyzIn.readline()).split(' ')
  # remove empty elements
  data = filter(None, data)
  aNameB.append(data[0])
  crd = [float(data[1]),float(data[2]),float(data[3])]
  aCrdB.append(crd)
#  Cb += array(crd)
aCrdB = vstack(aCrdB)
#Cb /= float(Nb)

# reading aligning atoms of B.xyz
itr += 1
while itr < len(sys.argv):
  listB.append(int(sys.argv[itr])-1)
  itr += 1
if((len(listB)<3) or (len(listB) != len(listA))):
  usage()


# reference atoms
rCrdA = aCrdA[listA][:]
rCrdB = aCrdB[listB][:]
Ca = mean(rCrdA,axis=0)
Cb = mean(rCrdB,axis=0)

# shift A.xyz and B.xyz to center
cCrdA = aCrdA - kron(Ca,ones((Na,1)))
cCrdB = aCrdB - kron(Cb,ones((Nb,1)))
crCrdA = rCrdA - kron(Ca,ones((len(listA),1)))
crCrdB = rCrdB - kron(Ca,ones((len(listB),1)))
# SVD and rotation
H = dot(transpose(crCrdA),crCrdB)
U, s, V = linalg.svd(H)
R = dot(transpose(V),transpose(U))
n1CrdA = transpose(dot(R,transpose(cCrdA)))+ kron(Cb,ones((Na,1)))

#print "A.xyz"
#print aCrdA
#print "B.xyz"
#print aCrdB
#print "rA"
#print rCrdA
#print "rB"
#print rCrdB
#print "crA"
#print crCrdA
#print "crB"
#print crCrdB

# print xyz format
printXYZ(Na,n1CrdA,aNameA)
#printXYZ(Nb,cCrdB,aNameB)
