#!/usr/bin/perl
use Cwd;
use Switch;
use strict;
use warnings;

# tuncation order 2<=C<=3
my $C=2;
# expansion order 4<=X<=8
my $X=4;
# spin dependent
#  S=1 for unpolarized
#  S=2 for polarized
my $S=1;

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

################################
# print CASINO corrlation.data #
# format for Jastrow factors   #
################################
# U terms and CHI terms are necessary
# while F terms could be omitted 
print <<EOF;
START JASTROW
Title
$ARGV[0]
Truncation order C
  $C
START U TERM
START SET 1
Spherical harmonic l,m
  0 0
Expansion order N_u
  $X
Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud)
  $S
Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)
  0                               0
Parameter values  ;  Optimizable (0=NO; 1=YES)
END SET 1
END U TERM
START CHI TERM
EOF

print "END CHI TERM\n";

###################################
# print atom type and coordinates #
###################################
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
