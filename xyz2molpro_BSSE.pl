#!/usr/bin/perl
use Cwd;
use Switch;
use List::Util qw(min max);
use POSIX;
use strict;
use warnings;

if(scalar(@ARGV) != 2){
  print "2 arguments needed: monomer1.xyz monomer2.xyz\n";
  die;
}

#my $PERIODIC_DISTANCE = $ARGV[1];
#my $PPDIR = '/home/chemie/changk/projects/03_Biotin/10_QMC/02_pwscfConvergenceTest/PP';
##my $PPDIR = './';
#my $CUTOFF = $ARGV[2];
#my $PPSTRING = '.pbe-hgh.UPF'."\n";
#my $FUNCTIONAL = 'PBE';

##################################
# read in number of atom $NA and #
# second line information $info2 #
##################################
open XYZ1, "<", "$ARGV[0]" or die $!;
open XYZ2, "<", "$ARGV[1]" or die $!;
my $NA1 = <XYZ1>;
my $NA2 = <XYZ2>;
$NA1 =~ s/\n//;
$NA2 =~ s/\n//;
my $NA3 = $NA1+$NA2;

#my $stem = $ARGV[0];
#$stem =~ s/\..*//;
my $info1 = <XYZ1>;
my $info2 = <XYZ2>;
#my $itr;
#
my ($Nve1, $Nve2, $Nve3) = (0,0,0);
my ($multi1, $multi2, $multi3) = (1,1,1);
my ($charge1, $charge2,$charge3) = (0,0,0);
my (@ATOMS1, @ATOMS2, @ATOMS3);
my (%GROUP);

##################################
# read atom type and coordinates #
##################################
for(my $itr=0;$itr<$NA1;$itr++){
  # read and string into componets
  my @input = split(' ',<XYZ1>);
  my @ATOM;
  my $AID = sprintf("%s%02d",$input[0],$itr+1);
  $ATOM[0] = atom_a2n($input[0]);
  for(my $i=1;$i<4;$i++){$ATOM[$i] = $input[$i];}
  # assign multiD data array
  push (@{$ATOMS1[$itr]}, @ATOM);
  $Nve1 += atom_ve($input[0]);
  $GROUP{$AID} = 1;
}
for(my $itr=0;$itr<$NA2;$itr++){
  # read and string into componets
  my @input = split(' ',<XYZ2>);
  my @ATOM;
  my $AID = sprintf("%s%02d",$input[0],$NA1+$itr+1);
  $ATOM[0] = atom_a2n($input[0]);
  for(my $i=1;$i<4;$i++){$ATOM[$i] = $input[$i];}
  # assign multiD data array
  push (@{$ATOMS2[$itr]}, @ATOM);
  $Nve2 += atom_ve($input[0]);
  $GROUP{$AID} = 2;
}
$Nve3 = $Nve1 + $Nve2;
$multi1 += $Nve1%2;
$multi2 += $Nve2%2;
$multi3 += $Nve3%2;
push(@ATOMS3,@ATOMS1,@ATOMS2);

if($multi1 != 1){
  $multi1 = 1;
  $charge1 = -1;
}
if($multi2 != 1){
  $multi2 = 1;
  $charge2 = -1;
}

#############################################
### read what every info at the end of file #
#############################################
##my @infoEND;
##for($itr=0;<XYZ>;$itr++){
##  #my @info = split(' ',$_);
##  push (@infoEND, $_);
##}
##my $Ni=$itr;
#
#
####################
## data processing #
####################
#my $pwd = cwd();
#
## numerical sort first column of @atoms
#@ATOMS1 = sort {$a->[0] <=> $b->[0]} @ATOMS1;
#@ATOMS2 = sort {$a->[0] <=> $b->[0]} @ATOMS2;
#
## atom type map
#my @aMAP;
#$aMAP[0]=0;
#my $NTYPE = 1;
#for(my $a=1;$a<$NA;$a++){
#  if($ATOMS[$a][0] != $ATOMS[$a-1][0]){
#    $aMAP[$NTYPE] = $a;
#    $NTYPE += 1;
#  }
#}
#push(@aMAP,$NA);
#
## shift to the corner of first quadrant
#for(my $i=1;$i<4;$i++){
#  my @vec = map { $_->[$i]} @ATOMS;
#  my $min = min @vec;
#  for(my $a=0;$a<$NA;$a++){
#    $ATOMS[$a][$i] -= $min;
#  }
#}
#
## shift molecule to MARGEIN from quadrant corner
#my @BOXD;
#for(my $i=1;$i<4;$i++){
#  my @vec = map { $_->[$i]} @ATOMS;
#  my $max = max @vec;
#  $BOXD[$i-1] = ceil($max + $PERIODIC_DISTANCE);
#}
#for(my $i=0;$i<3;$i++){
#  my $margin = $PERIODIC_DISTANCE/2.0;
#  for(my $a=0;$a<$NA;$a++){
#    $ATOMS[$a][$i+1] += $margin;
#  }
#}
#
#######################
# print Morpro header #
#######################
#memory,10000,m                   !memory in megawords, 1mw=7.6mb
#
print <<EOF;
memory,512,m                   !memory in megawords, 1mw=7.6mb

direct,dmp2=1
basis=6-31+G(D,P)

orient,noorient
angstrom
geometry={
EOF

###################################
# print atom type and coordinates #
###################################
# print XYZ for each atom type
for(my $i=0;$i<$NA3;$i++){
  my $name=sprintf("%s%02d",atom_a2n($ATOMS3[$i][0]),$i+1);
  print $name, " $ATOMS3[$i][1] $ATOMS3[$i][2] $ATOMS3[$i][3]\n";
}
print "}\n\n";
my ($dummy2,$dummy1);
for(my $i=0;$i<$NA1;$i++){
  my $name=sprintf("%s%02d",atom_a2n($ATOMS3[$i][0]),$i+1);
  $dummy1 .= ",$name";
}
for(my $i=$NA1;$i<$NA3;$i++){
  my $name=sprintf("%s%02d",atom_a2n($ATOMS3[$i][0]),$i+1);
  $dummy2 .= ",$name";
}
print <<EOF;
text, CP calculation for monomer1
dummy$dummy2; !monomer2 is now dummy
SET,CHARGE=$charge1
hf;           !scf for monomer1
{mp2;thresh,coeff=1d-2,energy=1d-4}
eny1=energy   !save mp2 energy as eny1

text, CP calculation for monomer2
dummy$dummy1; !monomer1 is now dummy
SET,CHARGE=$charge2
hf;           !scf for monomer2
{mp2;thresh,coeff=1d-2,energy=1d-4}
eny2=energy   !save mp2 energy as eny2

text, dimer calculation
dummy  !reset dummies
SET,CHARGE=$charge3
hf;           !scf for dimer
{mp2;thresh,coeff=1d-2,energy=1d-4}
eny3=energy   !save mp2 energy as eny3

text,computed interaction energies
de=(eny3-eny1-eny2)*tokcal    !CPC corrected interaction energy with fixed monomers
EOF

close(XYZ1);
close(XYZ2);

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

