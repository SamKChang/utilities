#!/usr/bin/perl
use Cwd;
use Switch;
use List::Util qw(min max);
use POSIX;
use strict;
use warnings;

if(scalar(@ARGV) != 4){
  print "4 arguments needed: xyz, margin, cutoff, ref.xyz\n";
  die;
}

my $PERIODIC_DISTANCE = 2*$ARGV[1];
my $margin = $ARGV[1];
#my $PPDIR = '/home/chemie/changk/projects/03_Biotin/10_QMC/02_pwscfConvergenceTest/PP';
#my $PPDIR = '/home/chemie/changk/packages/QuantumESPRESSO/PP/03_biotin/';
my $PPDIR = '/home/samio/Works/PhD/packages/QuantumESPRESSO/PP/03_biotin/';
#my $PPDIR = './';
my $CUTOFF = $ARGV[2];
my $PPSTRING = '_dcacp_blyp.UPF'."\n";
my $FUNCTIONAL = 'BLYP';

##################################
# read in number of atom $NA and #
# second line information $info2 #
##################################
open XYZ, "<", "$ARGV[0]" or die $!;
open XYZ1, "<", "$ARGV[3]" or die $!;
my $stem = $ARGV[0];
$stem =~ s/\..*//;

my $NA = <XYZ>;
$NA =~ s/\n//;
my $info2 = <XYZ>;

my $NA1 = <XYZ1>;
$NA1 =~ s/\n//;
my $info21 = <XYZ1>;

my $Nve = 0;
my $multi = 1;
my $charge = 0;
my @ATOMS;
my @ATOMS1;

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
#
# WARNING!!!!
# SPECIAL LOOP FOR BIOTIN PROJECT 
# WITH SPECIFIC NAMING SYSTEM
if(($ARGV[0] =~ m/rmv-..\.xyz/) and ($ARGV[0] !~ m/\-0[17]/)){
  $charge = -2;
}elsif($ARGV[0] =~ m/iso-07\.xyz/){
  $charge = -2;
}




############################
# read ref.xyz coordinates #
############################
for(my $itr=0;$itr<$NA1;$itr++){
  # read and string into componets
  my @input = split(' ',<XYZ1>);
  my @ATOM;
  $ATOM[0] = atom_a2n($input[0]);
  for(my $i=1;$i<4;$i++){$ATOM[$i] = $input[$i];}
  # assign multiD data array
  push (@{$ATOMS1[$itr]}, @ATOM);
}

###########################################
# read what every info at the end of file #
###########################################
my @infoEND;
my $Ni;
for($Ni=0;<XYZ>;$Ni++){
  #my @info = split(' ',$_);
  push (@infoEND, $_);
}


###################
# data processing #
###################
my $pwd = cwd();

# numerical sort first column of @atoms
@ATOMS = sort {$a->[0] <=> $b->[0]} @ATOMS;
@ATOMS1 = sort {$a->[0] <=> $b->[0]} @ATOMS1;

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
  my @vec = map { $_->[$i]} @ATOMS1;
  my $min = min @vec;
  for(my $a=0;$a<$NA;$a++){
    $ATOMS[$a][$i] -= $min;
  }
  for(my $a=0;$a<$NA1;$a++){
    $ATOMS1[$a][$i] -= $min;
  }
}

# shift molecule to MARGEIN from quadrant corner
my @BOXD;
for(my $i=1;$i<4;$i++){
  my @vec = map { $_->[$i]} @ATOMS1;
  my $max = max @vec;
  $BOXD[$i-1] = $max + $PERIODIC_DISTANCE;
}
for(my $i=0;$i<3;$i++){
  for(my $a=0;$a<$NA;$a++){
    $ATOMS[$a][$i+1] += $margin;
  }
}

# input_dft = '$FUNCTIONAL',
# ecutrho = 1000
######################
# print pwscf header #
######################
print <<EOF;
&control
 calculation = 'scf',
 prefix = '$stem',
 pseudo_dir = '$PPDIR',
 restart_mode='from_scratch',
 outdir='out',
 wf_collect=.true.
/
&system
 assume_isolated = 'mt',
 ibrav = 0,
 nat = $NA,
 ntyp = $NTYPE,
 ecutwfc = $CUTOFF,
 nosym = .true.,
EOF
if($charge!=0){
  print " tot_charge = $charge\n";
}
print <<EOF;
/
&electrons
 mixing_beta = 0.7,
 conv_thr = 1.D-9,
 diagonalization='davidson'
/
ATOMIC_SPECIES
EOF

###################################
# print atom type and coordinates #
###################################
# print XYZ for each atom type
for(my $t=0;$t<$NTYPE;$t++){
  # print pseudopotential files and atom names
  my $aType=atom_a2n($ATOMS[$aMAP[$t]][0]);
  my $mass;
  if($aType =~ m/^H/){
    $mass = 1;
  }else{
    $mass = $ATOMS[$aMAP[$t]][0] * 2;
  }
  print " $aType $mass ";
  print $aType,$PPSTRING;
}
#print "\nATOMIC_POSITIONS angstrom\n";
print "\nATOMIC_POSITIONS bohr\n";
for(my $a=0;$a<$NA;$a++){
  printf("  %-2s", atom_a2n($ATOMS[$a][0]));
  #for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS[$a][$i]);}
  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS[$a][$i]/0.52917721092);}
  print "\n";
}
print "K_POINTS gamma\n";
# print in bohr
print "CELL_PARAMETERS bohr\n";
printf(" %12.8f %12.8f %12.8f\n", $BOXD[0]/0.52917721092, 0, 0); 
printf(" %12.8f %12.8f %12.8f\n", 0, $BOXD[1]/0.52917721092, 0); 
printf(" %12.8f %12.8f %12.8f\n", 0, 0, $BOXD[2]/0.52917721092);
#print "CELL_PARAMETERS angstrom\n";
#printf(" %12.8f %12.8f %12.8f\n", $BOXD[0], 0, 0); 
#printf(" %12.8f %12.8f %12.8f\n", 0, $BOXD[1], 0); 
#printf(" %12.8f %12.8f %12.8f\n", 0, 0, $BOXD[2]);

#for(my $a=0;$a<$NA;$a++){
#  printf("%-2s", atom_a2n($ATOMS[$a][0]));
#  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS[$a][$i]);}
#  print "\n";
#}

#for(my $a=0;$a<$NA;$a++){
#  printf("%-2s", $ATOMS[$a][0]);
#  for(my $i=1;$i<4;$i++){printf(" % 11.8f",$ATOMS[$a][$i]);}
#  print "\n";
#}
#print "\n\n";
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

