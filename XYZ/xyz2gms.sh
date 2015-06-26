#!/bin/bash

cat $1|grep H|awk \
'{printf("  % .6f  % .6f  % .6f\n",$2,$3,$4)}' > H.crd
cat $1|grep C|awk \
'{printf("  % .6f  % .6f  % .6f\n",$2,$3,$4)}' > C.crd
cat $1|grep O|awk \
'{printf("  % .6f  % .6f  % .6f\n",$2,$3,$4)}' > O.crd
cat $1|grep N|awk \
'{printf("  % .6f  % .6f  % .6f\n",$2,$3,$4)}' > N.crd
cat $1|grep S|awk \
'{printf("  % .6f  % .6f  % .6f\n",$2,$3,$4)}' > S.crd
nH=$(cat H.crd|wc -l)
nC=$(cat C.crd|wc -l)
nO=$(cat O.crd|wc -l)
nN=$(cat N.crd|wc -l)
nS=$(cat S.crd|wc -l)

echo " \$CONTRL SCFTYP=ROHF RUNTYP=ENERGY MULT=1"
echo "    ISPHER=1 EXETYP=RUN MAXIT=200 ECP=READ DFTTYP=PBE \$END"
echo " \$SYSTEM MEMORY=150000000 \$END"
#echo " \$STATPT OPTTOL=1.0E-5  \$END"
#echo " \$BASIS  GBASIS=HW \$END"
echo " \$GUESS  GUESS=HUCKEL \$END"
echo " \$ECP"
if [ $nH -gt 0 ]; then
  cat H-qmc.ecp
  for i in `seq 2 $nH`;do
    echo "H-QMC"
  done
fi

if [ $nC -gt 0 ]; then
  cat C-qmc.ecp
  for i in `seq 2 $nC`;do
    echo "C-QMC"
  done
fi

if [ $nN -gt 0 ]; then
  cat N-qmc.ecp
  for i in `seq 2 $nN`;do
    echo "N-QMC"
  done
fi

if [ $nO -gt 0 ]; then
  cat O-qmc.ecp
  for i in `seq 2 $nO`;do
    echo "O-QMC"
  done
fi

if [ $nS -gt 0 ]; then
  cat S-qmc.ecp
  for i in `seq 2 $nS`;do
    echo "S-QMC"
  done
fi
echo " \$END"
echo " \$DATA"
echo $1
echo "C1"

if [ $nH -gt 0 ]; then
  for n in `seq 1 $(cat H.crd|wc -l)`;do 
    paste <(printf "H   1.0\n%.0s") <(sed -n "$n"p H.crd)
    cat H-qmc.basis
    echo ""
  done
fi

if [ $nC -gt 0 ]; then
  for n in `seq 1 $(cat C.crd|wc -l)`;do 
    paste <(printf "C   6.0\n%.0s") <(sed -n "$n"p C.crd)
    cat C-qmc.basis
    echo ""
  done
fi

if [ $nN -gt 0 ]; then
  for n in `seq 1 $(cat N.crd|wc -l)`;do 
    paste <(printf "N   7.0\n%.0s") <(sed -n "$n"p N.crd)
    cat N-qmc.basis
    echo ""
  done
fi

if [ $nO -gt 0 ]; then
  for n in `seq 1 $(cat O.crd|wc -l)`;do 
    paste <(printf "O   8.0\n%.0s") <(sed -n "$n"p O.crd)
    cat O-qmc.basis
    echo ""
  done
fi

if [ $nS -gt 0 ]; then
  for n in `seq 1 $(cat S.crd|wc -l)`;do 
    paste <(printf "S  16.0\n%.0s") <(sed -n "$n"p S.crd)
    cat S-qmc.basis
    echo ""
  done
fi

#if [ $nH -gt 0 ]; then
#  paste <(printf "H   1.0\n%.0s" $(seq 1 $nH)) H.crd
#  cat H-qmc.basis
#  echo ""
#fi
#
#if [ $nC -gt 0 ]; then
#  paste <(printf "C   6.0\n%.0s" $(seq 1 $nC)) C.crd
#  cat C-qmc.basis
#  echo ""
#fi
#
#if [ $nN -gt 0 ]; then
#  paste <(printf "N   7.0\n%.0s" $(seq 1 $nN)) N.crd
#  cat N-qmc.basis
#  echo ""
#fi
#
#if [ $nO -gt 0 ]; then
#  paste <(printf "O   8.0\n%.0s" $(seq 1 $nO)) O.crd
#  cat O-qmc.basis
#  echo ""
#fi
#
#if [ $nS -gt 0 ]; then
#  paste <(printf "S  16.0\n%.0s" $(seq 1 $nS)) S.crd
#  cat S-qmc.basis
#  echo ""
#fi

echo " \$END"

rm H.crd C.crd N.crd O.crd S.crd
