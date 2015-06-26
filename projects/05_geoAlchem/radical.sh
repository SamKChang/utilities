#!/bin/bash

N0=`sed -n 1p $1`
tmp=`sed -n 2p $1`
N1=`echo $tmp-1|bc`
N2=`echo "$N0 - $N1"|bc`

L1=`echo $tmp+2|bc -l`
D2=`echo $L1+1|bc`
L2=`cat $1|wc -l`
S2=`seq -s " " $D2 $L2`
Name1="ra-"$1
Name2="rb-"$1

echo $N1 > $Name1
echo "" >> $Name1
sed -n "4,$L1 p" $1 >> $Name1
echo $N2 > $Name2
echo "" >> $Name2
sed -n 3p $1 >> $Name2
sed -n "$D2,$L2 p" $1 >> $Name2

#echo $N1 # > $Name1
#echo "" # >> $Name1
#sed -n "4,$L1 p" $1 # >> $Name1
#echo $N2 # > $Name2
#echo "" # >> $Name2
#sed -n 3p $1 # >> $Name2
#sed -n "$D2,$L2 p" $1 # >> $Name2
