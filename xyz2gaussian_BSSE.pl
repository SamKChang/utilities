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

my $stem = $ARGV[0];
$stem =~ s/[AB]\..*//;
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
#$multi3 += $Nve3%2;
push(@ATOMS3,@ATOMS1,@ATOMS2);

if($ARGV[0] =~ m/^an/){
  if($multi1 != 1){
    $multi1 = 1;
    $charge1 = -1;
  }
  if($multi2 != 1){
    $multi2 = 1;
    $charge2 = -1;
  }
}else{
  if($multi1 != 1){
    $multi1 = 1;
    $charge1 = 1;
  }
  if($multi2 != 1){
    $multi2 = 1;
    $charge2 = 1;
  }
}
$charge3 = $charge1 + $charge2;

###########################
# print Gaussian09 header #
###########################
print <<EOF;
%nproc=1
# pbepbe/aug-cc-pVTZ Counterpoise=2 EmpiricalDispersion=GD3BJ 6d 10f nosymm Scf(maxcycle=1000,vshift,verytight) int(grid=ultrafine)

$stem

$charge3   $multi3   $charge1   $multi1   $charge2   $multi2
EOF

###################################
# print atom type and coordinates #
###################################
# print XYZ for each atom type
for(my $i=0;$i<$NA1;$i++){
  my $name=sprintf("%s(Fragment=1)",atom_a2n($ATOMS1[$i][0]),$i+1);
  print "$name ";
  for(my $j=1;$j<4;$j++){
    printf("% 12.8f ", $ATOMS1[$i][$j]);
  }
  print "\n";
}
for(my $i=0;$i<$NA2;$i++){
  my $name=sprintf("%s(Fragment=2)",atom_a2n($ATOMS2[$i][0]),$i+1);
  print "$name ";
  for(my $j=1;$j<4;$j++){
    printf("% 12.8f ", $ATOMS2[$i][$j]);
  }
  print "\n";
}
print "\n\n";
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

