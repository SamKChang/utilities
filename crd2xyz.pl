#!/usr/bin/perl
use Cwd;
use Switch;
use strict;
use warnings;

##################################
# read in number of atom $NA and #
# second line information $info2 #
##################################
open CRD, "<", "$ARGV[0]" or die $!;
my $stem = $ARGV[0];
$stem =~ s/\..*//;
my @ATOMS = <CRD>;
my $itr=0;
my @XYZ;
my $A=0.52917721092;

foreach(@ATOMS){
  my @data = split(" ",$ATOMS[$itr]);
  my $str= sprintf("%-2s % 10.8f % 10.8f % 10.8f\n", atom_a2n(int($data[0])), $data[1]*$A, $data[2]*$A,$data[3]*$A);
  push(@XYZ,$str);
  $itr +=1 ;
}

print $itr,"\n\n";
print @XYZ;


###################################
## read atom type and coordinates #
###################################
#for(my $itr=0;$itr<$NA;$itr++){
#  # read and string into componets
#  my @ATOM = split(' ',<XYZ>);
#  # assign multiD data array
#  push (@{$ATOMS[$itr]}, @ATOM);
#  $Nve += atom_ve($ATOM[0]);
#}
#$multi += $Nve%2;
#
## set multiplicity to 1 by adding electrons
#if($multi != 1){
#  $multi = 1;
#  $charge = -1;
#}
#
############################################
## read what every info at the end of file #
############################################
#my @infoEND;
#for($itr=0;<XYZ>;$itr++){
#  #my @info = split(' ',$_);
#  push (@infoEND, $_);
#}
#my $Ni=$itr;
#
#my $pwd = cwd();
#
#
#####################
## Gaussian09 NOTES #
#####################
## command for Hamiltonian components print out
## IOp(3/5=17) request MWB default basis
## IOp(3/6=6) request MWB default ECP
## IOp(3/33=4) request print out 1-particle Matrix components
## IOp(2/12=3) supress "atoms too near" check
## guess(indo) request semiemperical initial guess for density matrix
#
########################
## Gaussian09 commands #
########################
## hf 6d 10f nosymm Scf(maxcycle=1,vshift,verytight) int(grid=ultrafine) IOp(3/5=17) IOp(3/6=6) IOp(3/33=4) IOp(2/12=3) guess(indo)
## ropbe1pbe 6d 10f nosymm Scf(maxcycle=1000,vshift,verytight) int(grid=ultrafine) IOp(3/5=17) IOp(3/6=6) IOp(2/12=3)
## ropbe1pbe 6d 10f nosymm Scf(maxcycle=1000,vshift,verytight) int(grid=ultrafine) IOp(3/5=17) IOp(3/6=6) IOp(3/33=4) guess(indo)
##%chk=$pwd/com/$stem
#
############################
## print Gaussian09 header #
############################
#print <<EOF;
#%nproc=1
## mp2/aug-cc-pVTZ 6d 10f nosymm Scf(maxcycle=1000,vshift,verytight) int(grid=ultrafine)
#
#$stem
#
#$charge   $multi
#EOF
#
####################################
## print atom type and coordinates #
####################################
#for(my $a=0;$a<$NA;$a++){
#  printf("%-2s", $ATOMS[$a][0]);
#  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS[$a][$i]);}
#  print "\n";
#}
#print "\n\n";
#close(XYZ);
#
####################
## useful routines #
###################
## count valence electrons
#sub atom_ve{
#  switch($_[0]){
#    case "H" {return 1}
#    case "C" {return 4}
#    case "Si"{return 4}
#    case "Ge"{return 4}
#    case "Sn"{return 4}
#    case "N" {return 5}
#    case "P" {return 5}
#    case "As"{return 5}
#    case "Sb"{return 5}
#    case "O" {return 6}
#    case "S" {return 6}
#    case "Se"{return 6}
#    case "Te"{return 6}
#    case "F" {return 7}
#    case "Cl"{return 7}
#    case "Br"{return 7}
#    case "I" {return 7}
#    else     {
#      print "element: $_[0] not found\n";
#      die;
#    }
#  }
#}
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
