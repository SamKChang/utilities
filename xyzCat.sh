#!/bin/bash

N1=$(head -n 1 $1)
N2=$(head -n 1 $2)
N=$(($N1 + $N2))

echo $N
echo ""
tail -n $N1 $1
tail -n $N2 $2
