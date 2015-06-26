#!/bin/bash

cat $1|grep H|awk \
'{printf("  % .14f  % .14f  % .14f\n",$2,$3,$4)}' > H.crd
cat $1|grep C|awk \
'{printf("  % .14f  % .14f  % .14f\n",$2,$3,$4)}' > C.crd
cat $1|grep O|awk \
'{printf("  % .14f  % .14f  % .14f\n",$2,$3,$4)}' > O.crd
cat $1|grep N|awk \
'{printf("  % .14f  % .14f  % .14f\n",$2,$3,$4)}' > N.crd

nH=$(cat H.crd|wc -l)
nC=$(cat C.crd|wc -l)
nO=$(cat O.crd|wc -l)
nN=$(cat N.crd|wc -l)

echo "&CPMD"
echo " OPTIMIZE WAVEFUNCTION"
echo " BENCHMARK"
echo "  1 0 0 0 0 0 0 0 0 0"
echo " MIRROR"
echo " MAXITER"
echo "  1.0e5"
#echo " OPTIMIZER PCG"
echo "&END"
echo ""
echo "&DFT"
echo " FUNCTIONAL BLYP"
echo "&END"
echo ""
echo "&SYSTEM"
echo " SYMMETRY"
echo "  ISOLATED"
echo " POISSON SOLVER TUCKERMAN"
echo " ANGSTROM"
echo " CELL ABSOLUTE"
echo "  22 22 22 0 0 0"
echo " CUTOFF"
echo "  110"
echo " MESH"
echo "  220 220 220"
echo " KPOINTS MONKHORST-PACK"
echo "  1 1 1"
if [ -e ".CHRG" ]; then
	chrg=$(cat .CHRG)
	echo " CHARGE"
	echo "  $chrg"
fi
echo "&END"
echo ""
echo "&ATOMS"

if [ $nH -gt 0 ]; then
#	echo "*H_blyp_mt.psp KLEINMAN-BYLANDER"
#	echo " LMAX=P"
	echo "*H_dcacp_blyp_mt.psp KLEINMAN-BYLANDER"
	echo " LMAX=F LOC=P SKIP=D"
#	echo "*H_dcacp_blyp_gs.psp"
#	echo " LMAX=F"
	echo "  $nH"
	cat H.crd
	echo ""
fi

if [ $nC -gt 0 ]; then
#	echo "*C_blyp_mt.psp KLEINMAN-BYLANDER"
#	echo " LMAX=D"
	echo "*C_dcacp_blyp_mt.psp KLEINMAN-BYLANDER"
	echo " LMAX=F LOC=D"
#	echo "*C_dcacp_blyp_gs.psp"
#	echo " LMAX=F"
	echo "  $nC"
	cat C.crd
	echo ""
fi

if [ $nN -gt 0 ]; then
#	echo "*N_blyp_mt.psp KLEINMAN-BYLANDER"
#	echo " LMAX=D"
	echo "*N_dcacp_blyp_mt.psp KLEINMAN-BYLANDER"
	echo " LMAX=F LOC=D"
#	echo "*N_dcacp_blyp_gs.psp"
#	echo " LMAX=F"
	echo "  $nN"
	cat N.crd
	echo ""
fi

if [ $nO -gt 0 ]; then
#	echo "*O_blyp_mt.psp KLEINMAN-BYLANDER"
#	echo " LMAX=D"
	echo "*O_dcacp_blyp_mt.psp KLEINMAN-BYLANDER"
	echo " LMAX=F LOC=P SKIP=D"
#	echo "*O_dcacp_blyp_gs.psp"
#	echo " LMAX=F"
	echo "  $nO"
	cat O.crd
	echo ""
fi
echo "&END"

rm H.crd C.crd N.crd O.crd
