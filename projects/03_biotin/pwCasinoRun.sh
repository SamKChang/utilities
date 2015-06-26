#!/bin/bash

# pwscf/casino calculations from xyz structures
# expect structrue folders as input argument
# following perl script for input interfaces:
#   xyz2pw_center.pl XYZ boxSize(Ang) cutoff (Ry)
#   xyz2casino-corr.pl XYZ
#   xyz2casino-Jopt.pl XYZ VMCStep_opt
#   xyz2casino-dmc.pl XYZ VMCStep DMCStep DMCdt

# system setup
NCPU=$2
PWEXE="pw.x"
RPATH=$PWD
XYZPATH=`readlink -f $1`
PP_PATH="/home/samio/Works/PhD/packages/CASINO/PP/"
MPISTR="mpirun -np $NCPU"
#MPISTR="mpirun -np $NCPU -mca btl tcp,self"

# computation setup
boxSize="21"
cutoff="250"
VMCStep_opt="10000"
VMCStep="5000"
DMCStep="100000"
DMCdt=(0.010 0.012 0.014 0.016)

# loop through all xyz files in target folder
for f in `ls $XYZPATH|grep "\.xyz$"|head -n 1`;do
#for f in `ls $XYZPATH|grep "anDIR_26-00.*\.xyz$"`;do
  FILE=${f%.*}
  INP=${f%.*}.inp
  mkdir $FILE
  cd $FILE

    # PWSCF for blip wavefunction
    mkdir pwscf
    cd pwscf
      cp $XYZPATH/$f .
      xyz2pw_center.pl $f $boxSize $cutoff > $INP
      mkdir out
      cd out
        echo "&inputpp"                 >> pw2casino.dat
#        echo "blip_multiplicity=2.d0"   >> pw2casino.dat
        echo "blip_single_prec=.true."  >> pw2casino.dat
        echo "/                      "  >> pw2casino.dat
      cd ..
      $MPISTR $PWEXE -pw2casino < $INP > ${f%.*}.out
    cd ..

    # Jastrow factor optimization
    mkdir vmc_opt
    cd vmc_opt
      cp $XYZPATH/$f .
      cp $PP_PATH/*.data .
      ln $RPATH/$FILE/pwscf/out/*.bwfn.data* bwfn.data.b1
      xyz2casino-Jopt.pl $f $VMCStep_opt > input
      xyz2casino-corr.pl $f > correlation.data
      runqmc -p $NCPU
      EMIN=`ve|awk '/VMC #/{print $5}'\
            |sed 's/(.*//g'|sort -n|tail -n 1\
            |grep -oh "[-0-9\.]*"|tail -n 1|sed 's/-/\\\-/g'`
      VMIN=`ve|grep "VMC #"|grep $EMIN|awk '{print $9}'\
            |sed 's/(.*)//g'|sort -n|head -n 1\
            |grep -oh "[-0-9\.]*"|tail -n 1`
      OPT=`ve|grep "VMC #.*$EMIN.*$VMIN"\
           |awk '{print $10}'|sed 's/[()]//g'`
      rm bwfn.data.b1
    cd ..

    # DMC statistic collection with zero time-step extrapolation
    #(start from scratch with optimized J-factors)
    for dt in ${DMCdt[*]};do
      dmcdt=`echo $dt|sed 's/\.//g'`
      dmc="dmc$dmcdt"
      mkdir $dmc
      cd $dmc
        cp $XYZPATH/$f .
        cp $PP_PATH/*.data .
        ln $RPATH/$FILE/pwscf/out/*.bwfn.data* bwfn.data.b1
        cp $RPATH/$FILE/vmc_opt/$OPT correlation.data
        xyz2casino-dmc.pl $f $VMCStep $DMCStep $dt > input
        runqmc -p $NCPU
        rm bwfn.data.b1
      cd ..
    done
  cd ..
  rm $RPATH/$FILE/pwscf/out/*.bwfn.data*
done
