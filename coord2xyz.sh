#!/bin/bash

N=$(cat $1|wc -l)
n=$(echo "$N - 3"|bc)
n1=$(echo "$n + 1"|bc)

head -n $n1 $1 | tail -n $n \
|awk '{printf("% .12f % .12f % .12f\n",$1*0.529,$2*0.529,$3*0.529)}'\
>xyz

head -n $n1 $1 | tail -n $n \
|awk '{print $4}' > atom

cat atom | tr '[a-z]' '[A-Z]' > ATOM

echo -e "$n\n"
paste ATOM xyz

rm xyz atom ATOM
