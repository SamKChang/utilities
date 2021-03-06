#!/usr/bin/perl
use Cwd;
use Switch;
use List::Util qw(min max);
use POSIX;
use strict;
use warnings;

# tuncation order 2<=C<=3
my $C=3;
# expansion order 4<=X<=8
my $X=6;
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

#### shift to the corner of first quadrant
###for(my $i=1;$i<4;$i++){
###  my @vec = map { $_->[$i]} @ATOMS;
###  my $min = min @vec;
###  for(my $a=0;$a<$NA;$a++){
###    $ATOMS[$a][$i] -= $min;
###  }
###}
###
#### shift molecule to MARGEIN from quadrant corner
###my @BOXD;
###for(my $i=1;$i<4;$i++){
###  my @vec = map { $_->[$i]} @ATOMS;
###  my $max = max @vec;
###  $BOXD[$i-1] = $max + $PERIODIC_DISTANCE;
###}
###for(my $i=0;$i<3;$i++){
###  for(my $a=0;$a<$NA;$a++){
###    $ATOMS[$a][$i+1] += $margin;
###  }
###}
###


#########################
#### print pwscf header #
#########################
###print <<EOF;
###&control
### calculation = 'scf',
### prefix = '$stem',
### pseudo_dir = '$PPDIR',
### restart_mode='from_scratch',
### outdir='out',
### wf_collect=.false.
###/
###&system
### assume_isolated = 'mt',
### ibrav = 0,
### nat = $NA,
### ntyp = $NTYPE,
### ecutwfc = $CUTOFF,
### nosym = .true.,
### noinv =.true.,
###EOF
#### use default dft functional from PPs
#### input_dft = '$FUNCTIONAL',
####  print " tot_charge = -2,\n";
###if($charge!=0){
###  print " tot_charge = $charge,\n";
###}
###print <<EOF;
###/
###&electrons
### mixing_beta = 0.7, diagonalization = 'david', conv_thr =  1.0d-9,
###/
###ATOMIC_SPECIES
###EOF

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
Number of sets
 1
START SET 1
Spherical harmonic l,m
 0 0
Expansion order N_u
 $X
Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud)
 $S
Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)
 0                               1
Parameter values  ;  Optimizable (0=NO; 1=YES)
END SET 1
END U TERM
START CHI TERM
Number of sets
 $NTYPE
EOF
for(my $t=0;$t<$NTYPE;$t++){
  my $ti = $t+1;
  print <<EOF; 
START SET $ti
Spherical harmonic l,m
 0 0
Number of atoms in set
EOF
  print " ",$aMAP[$t+1]-$aMAP[$t],"\n";
  print "Labels of the atoms in this set\n";
  for(my $i=$aMAP[$t];$i<$aMAP[$t+1];$i++){print " ", $i+1}
  print "\n";
  print <<EOF;
Impose electron-nucleus cusp (0=NO; 1=YES)
 0
Expansion order N_chi
 $X
Spin dep (0->u=d; 1->u/=d)
 $S
Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)
 0                               1
Parameter values  ;  Optimizable (0=NO; 1=YES)
END SET $ti
EOF
}
print <<EOF;
END CHI TERM
START F TERM
Number of sets
 $NTYPE
EOF
for(my $t=0;$t<$NTYPE;$t++){
  my $ti = $t+1;
  print <<EOF; 
START SET $ti
Number of atoms in set
EOF
  print " ",$aMAP[$t+1]-$aMAP[$t],"\n";
  print "Labels of the atoms in this set\n";
  for(my $i=$aMAP[$t];$i<$aMAP[$t+1];$i++){print " ", $i+1}
  print "\n";
  print <<EOF;
Prevent duplication of u term (0=NO; 1=YES)
 0
Prevent duplication of chi term (0=NO; 1=YES)
 0
Electron-nucleus expansion order N_f_eN
 2
Electron-electron expansion order N_f_ee
 2
Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud)
 $S
Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)
 0                   1
Parameter values  ;  Optimizable (0=NO; 1=YES)
END SET $ti
EOF
}
print <<EOF;
END F TERM
END JASTROW
EOF
######################################
#### print atom type and coordinates #
######################################
#### print XYZ for each atom type
###for(my $t=0;$t<$NTYPE;$t++){
###  # print pseudopotential files and atom names
###  my $aType=atom_a2n($ATOMS[$aMAP[$t]][0]);
###  my $mass;
###  if($aType =~ m/^H/){
###    $mass = 1;
###  }else{
###    $mass = $ATOMS[$aMAP[$t]][0] * 2;
###  }
###  print " $aType $mass ";
###  print $aType,$PPSTRING;
###}
###print "\nATOMIC_POSITIONS bohr\n";
####print "\nATOMIC_POSITIONS angstrom\n";
###for(my $a=0;$a<$NA;$a++){
###  printf("  %-2s", atom_a2n($ATOMS[$a][0]));
###  # print Bohr
###  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS[$a][$i]/0.52917721092);}
###  # print Angstrom
###  #for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS[$a][$i]);}
###  print "\n";
###}
###print "K_POINTS gamma\n";
###print "CELL_PARAMETERS bohr\n";
###printf(" %12.8f %12.8f %12.8f\n", $BOXD[0]/0.52917721092, 0, 0);
###printf(" %12.8f %12.8f %12.8f\n", 0, $BOXD[1]/0.52917721092, 0);
###printf(" %12.8f %12.8f %12.8f\n", 0, 0, $BOXD[2]/0.52917721092);
####print "CELL_PARAMETERS angstrom\n";
####printf(" %12.8f %12.8f %12.8f\n", $BOXD[0], 0, 0);
####printf(" %12.8f %12.8f %12.8f\n", 0, $BOXD[1], 0);
####printf(" %12.8f %12.8f %12.8f\n", 0, 0, $BOXD[2]);
###
####for(my $a=0;$a<$NA;$a++){
####  printf("%-2s", atom_a2n($ATOMS[$a][0]));
####  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS[$a][$i]);}
####  print "\n";
####}
###
####for(my $a=0;$a<$NA;$a++){
####  printf("%-2s", $ATOMS[$a][0]);
####  for(my $i=1;$i<4;$i++){printf(" % 11.8f",$ATOMS[$a][$i]);}
####  print "\n";
####}
####print "\n\n";
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

