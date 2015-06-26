#!/usr/bin/perl
use Cwd;
use Switch;
use List::Util qw(min max);
use POSIX;
use strict;
use warnings;

my $lambda = $ARGV[2];

my $margin = $ARGV[3];
my $PERIODIC_DISTANCE = 2*$margin;
my $CUTOFF = 50;
my $FUNCTIONAL = 'PBE';
#my $ppFUNCTION = lc $FUNCTIONAL . "_dcacp";
#my $ppFUNCTION = lc $FUNCTIONAL;
my $PPString = "MT_BLYP_DCACP.psp KLEINMAN-BYLANDER";

##################################
# read in number of atom $NA and #
# second line information $info2 #
##################################
open XYZ, "<", "$ARGV[0]" or die $!;
my $stem = $ARGV[0];
$stem =~ s/\..*//;
my $NA = <XYZ>;
$NA =~ s/\n//;
my $info2 = <XYZ>;
my $Nve = 0;
my $multi = 1;
my $charge = 0;
my @ATOMS;
my $itr;

##################################
# read atom type and coordinates #
##################################
for(my $itr=0;$itr<$NA;$itr++){
  # read and string into componets
  my @input = split(' ',<XYZ>);
  my @ATOM;
  $ATOM[0] = atom_a2n($input[0]);
  for(my $i=1;$i<4;$i++){$ATOM[$i] = $input[$i];}
  # assign multiD data array
  push (@{$ATOMS[$itr]}, @ATOM);
  $Nve += atom_ve($input[0]);
}
$multi += $Nve%2;
if($multi != 1){ 
  $charge = -1; 
  $multi = 1;
}

###################
# second XYZ file #
###################
open XYZ2, "<", "$ARGV[1]" or die $!;
my $stem2 = $ARGV[1];
$stem2 =~ s/.*\/(.*)\..*/$1/;
$stem2 =~ s/\..*//;
my $NA2 = <XYZ2>;
my $info22 = <XYZ2>;
my $Nve2 = 0;
my $multi2 = 1;
my $charge2 = 0;
my @ATOMS2;

for(my $itr=0;$itr<$NA2;$itr++){
  # read and string into componets
  my @ATOM2 = split(' ',<XYZ2>);
  # assign multiD data array
  push (@{$ATOMS2[$itr]}, @ATOM2);
  $Nve2 += atom_ve($ATOM2[0]);
}
$multi2 += $Nve2%2;

# set multiplicity to 1 by adding electrons
if($multi2 != 1){
  $multi2 = 1;
  $charge2 = -1;
}



###########################################
# read what every info at the end of file #
###########################################
my @infoEND;
for($itr=0;<XYZ>;$itr++){
  #my @info = split(' ',$_);
  push (@infoEND, $_);
}
my $Ni=$itr;


###################
# data processing #
###################
my $pwd = cwd();

# numerical sort first column of @atoms
@ATOMS = sort {$a->[0] <=> $b->[0]} @ATOMS;

# atom type map
my @aMAP;
$aMAP[0]=0;
my $NTYPE = 1;
for(my $a=1;$a<$NA;$a++){
  if($ATOMS[$a][0] != $ATOMS[$a-1][0]){
    $aMAP[$NTYPE] = $a;
    $NTYPE += 1;
  }
}
push(@aMAP,$NA);

# shift to the corner of first quadrant
for(my $i=1;$i<4;$i++){
  my @vec = map { $_->[$i]} @ATOMS;
  my $min = min @vec;
  for(my $a=0;$a<$NA;$a++){
    $ATOMS[$a][$i] -= $min;
  }
}

# shift molecule to MARGEIN from quadrant corner
my @BOXD;
for(my $i=1;$i<4;$i++){
  my @vec = map { $_->[$i]} @ATOMS;
  my $max = max @vec;
  $BOXD[$i-1] = $max + $PERIODIC_DISTANCE;
}
for(my $i=0;$i<3;$i++){
  for(my $a=0;$a<$NA;$a++){
    $ATOMS[$a][$i+1] += $margin;
  }
}

#####################
# print cpmd header #
#####################
#-------------#
# initial run #
#-------------#
# OPTIMIZE WAVEFUNCTION
# CONVERGENCE ORBITALS
#  1.0E-7
# MAXITER
#  1E5
# CENTER MOLECULE OFF
# MIRROR

#-----------------------#
# KS orbital generation #
#-----------------------#
# RESTART WAVEFUNCTION
# KOHN-SHAM ENERGY
#  185
# LANCZOS PARAMETER N=5
#  1000  16  20  1.D-9
#  0.05          1.D-11
#  0.01          1.D-13
#  0.0025        1.D-16
#  0.001         1.D-18
# CENTER MOLECULE OFF
# RHOOUT
# MIRROR
# MAXITER
#  1E5

#----------------#
# 1st order scan #
#----------------#
# RESTART WAVEFUNCTION
# OPTIMIZE WAVEFUNCTION
# MIRROR
# MAXITER
#  1
# CENTER MOLECULE OFF
# BENCHMARK
#  1 0 0 0 0 0 0 0 0 0

#----------------#
# 2nd order scan #
#----------------#
# RESTART WAVEFUNCTION
# MIRROR
# MAXITER
#  1
# CENTER MOLECULE OFF
# BENCHMARK
#  1 0 0 0 0 0 0 0 0 0
# KOHN-SHAM ENERGIES
#  185
# LANCZOS PARAMETER
#  0  16  20  1

# KANCZOS DIAGONALIZATION RESET=2
# CONVERGENCE ORBITALS
#  1.0E-7
# BENCHMARK
#  1 0 0 0 0 0 0 0 0 0

# SYMMETRY
#  ORTHORHOMBIC

# OPTIMIZE GEOMETRY
print <<EOF;
&INFO
 box margin: $margin
 ref: $ARGV[0]
 tar: $ARGV[1]
 lambda: $lambda
&END

&CPMD
 OPTIMIZE GEOMETRY
 CONVERGENCE ORBITALS
  1.0E-4
 MAXITER
  1E5
 CENTER MOLECULE OFF
 MIRROR
 RHOOUT
&END

&DFT
 FUNCTIONAL $FUNCTIONAL
&END

&SYSTEM
 SYMMETRY
  ISOLATED
 POISSON SOLVER TUCKERMAN
 CELL ABSOLUTE
  $BOXD[0]  $BOXD[1]  $BOXD[2] 0 0 0
 ANGSTROM
 CUTOFF
  $CUTOFF
EOF
if($charge!=0){
  print " CHARGE\n  $charge\n";
}
print <<EOF;
&END

&ATOMS
EOF


###################################
# print atom type and coordinates #
###################################
for(my $t=0;$t<$NTYPE;$t++){
  my $NAME_A = atom_a2n($ATOMS[$aMAP[$t]][0]);
  my $NT = $aMAP[$t+1] - $aMAP[$t];
  #print "*$NAME_A","_q",atom_ve("$NAME_A"),"_$ppFUNCTION",".psp ","FRAC\n";
  print "*$NAME_A",,"_$PPString"," FRAC\n";
  print " LMAX=F\n";
  print "  $NT\n";
  for(my $a=$aMAP[$t];$a<$aMAP[$t+1];$a++){
    print "  ";
    for(my $i=1;$i<4;$i++){
      printf("% 12.8f ",$ATOMS[$a][$i]);
    }
    print "\n";
  }
  print "\n";
}

####################
# print constraint #
####################
#print <<EOF;
#CONSTRAINTS
# FIX ELEMENT
#  6
# FIX ELEMENT
#  7
# FIX ELEMENT
#  8
# FIX ELEMENT
#  16
#END CONSTRAINTS
#EOF


print "&END\n";
close(XYZ);

###################
# useful routines #
##################
# count valence electrons
sub atom_a2n{
  switch($_[0]){
    case "H" {return  1}
    case "B" {return  5}
    case "C" {return  6}
    case "N" {return  7}
    case "O" {return  8}
    case "F" {return  9}
    case "Si"{return 14}
    case "P" {return 15}
    case "S" {return 16}
    case "Cl"{return 17}
    case "Ge"{return 32}
    case "As"{return 33}
    case "Se"{return 34}
    case "Br"{return 35}
    case "Sn"{return 50}
    case "Sb"{return 51}
    case "Te"{return 52}
    case "I" {return 53}
    case "1" {return "H" }
    case "5" {return "B" }
    case "6" {return "C" }
    case "7" {return "N" }
    case "8" {return "O" }
    case "9" {return "F" }
    case "14" {return "Si"}
    case "15" {return "P" }
    case "16" {return "S" }
    case "17" {return "Cl"}
    case "32" {return "Ge"}
    case "33" {return "As"}
    case "34" {return "Se"}
    case "35" {return "Br"}
    case "50" {return "Sn"}
    case "51" {return "Sb"}
    case "52" {return "Te"}
    case "53" {return "I" }
    #else    {
    #  print "element: $_[0] not found\n";
    #  die;
    #}
  }
}

sub atom_ve{
  switch($_[0]){
    case "H" {return 1}
    case "B" {return 3}
    case "C" {return 4}
    case "Si"{return 4}
    case "Ge"{return 4}
    case "Sn"{return 4}
    case "N" {return 5}
    case "P" {return 5}
    case "As"{return 5}
    case "Sb"{return 5}
    case "O" {return 6}
    case "S" {return 6}
    case "Se"{return 6}
    case "Te"{return 6}
    case "F" {return 7}
    case "Cl"{return 7}
    case "Br"{return 7}
    case "I" {return 7}
    else     {
      print "element: $_[0] not found\n";
      die;
    }
  }
}
