#!/usr/bin/perl
use Cwd;
use Switch;
use strict;
use warnings;

##################################
# read in number of atom $NA and #
# second line information $info2 #
##################################
open XYZ, "<", "$ARGV[0]" or die $!;
my $stem = $ARGV[0];
$stem =~ s/\..*//;
my $NA = <XYZ>;
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
  my @ATOM = split(' ',<XYZ>);
  # assign multiD data array
  push (@{$ATOMS[$itr]}, @ATOM);
  $Nve += atom_ve($ATOM[0]);
}
$multi += $Nve%2;

# set multiplicity to 1 by adding electrons
if($multi != 1){
  $multi = 1;
  $charge = -1;
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

my $pwd = cwd();


####################
# Gaussian09 NOTES #
####################
# command for Hamiltonian components print out
# IOp(3/5=17) request MWB default basis
# IOp(3/6=6) request MWB default ECP
# IOp(3/33=4) request print out 1-particle Matrix components
# IOp(2/12=3) supress "atoms too near" check
# guess(indo) request semiemperical initial guess for density matrix

#######################
# Gaussian09 commands #
#######################
# hf 6d 10f nosymm Scf(maxcycle=1,vshift,verytight) int(grid=ultrafine) IOp(3/5=17) IOp(3/6=6) IOp(3/33=4) IOp(2/12=3) guess(indo)
# ropbe1pbe 6d 10f nosymm Scf(maxcycle=1000,vshift,verytight) int(grid=ultrafine) IOp(3/5=17) IOp(3/6=6) IOp(2/12=3)
# ropbe1pbe 6d 10f nosymm Scf(maxcycle=1000,vshift,verytight) int(grid=ultrafine) IOp(3/5=17) IOp(3/6=6) IOp(3/33=4) guess(indo)
#%chk=$pwd/com/$stem
# ropbe1pbe 6d 10f nosymm Scf(maxcycle=1000,vshift,verytight) int(grid=superfine) IOp(3/5=17) IOp(3/6=6) IOp(3/33=4) guess(indo)

###########################
# print Gaussian09 header #
###########################
print <<EOF;
%nproc=1
%chk=$pwd/$stem
# pbepbe GEN 6d 10f nosymm Scf(maxcycle=1000,verytight) int(grid=ultrafine)

$stem

$charge   $multi
EOF

###################################
# print atom type and coordinates #
###################################
for(my $a=0;$a<$NA;$a++){
  printf("%-2s", $ATOMS[$a][0]);
  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS[$a][$i]);}
  print "\n";
}
print "\n\n";
close(XYZ);

###################
# useful routines #
##################
# count valence electrons
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
