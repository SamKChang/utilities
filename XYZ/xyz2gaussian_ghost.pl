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
$stem =~ s/.*\/(.*)\..*/$1/;
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
# pbepbe/def2-TZVP 6d 10f nosymm Scf(maxcycle=1000,verytight) int(grid=ultrafine)

###########################
# print Gaussian09 header #
###########################
print <<EOF;
%nproc=1
%chk=$pwd/$stem
# pbepbe gen 6d 10f nosymm Scf(maxcycle=1000,verytight) int(grid=ultrafine)

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
for(my $d=0.5;$d<=3.1;$d+=0.1){
  if(abs($d-$ATOMS[0][1])>0.01){
    printf("%-2s", "H-Bq");
    printf(" % 12.8f",$d);
    printf(" % 12.8f",0);
    printf(" % 12.8f",0);
    print "\n";
  }
}

my $AStr = $ATOMS[1][0];
print <<EOF;

$AStr     0
S   1   1.00
  45069.4640220              1.0000000
S   1   1.00                          
   6755.9768656              1.0000000
S   1   1.00                          
   1537.6502864              1.0000000
S   1   1.00                          
    435.51697667             1.0000000
S   1   1.00                          
    142.28655638             1.0000000
S   1   1.00                          
     51.692153804            1.0000000
S   1   1.00                          
     20.315870490            1.0000000
S   1   1.00                          
      8.2021942646           1.0000000
S   1   1.00
      1.9681276278           1.0000000        
S   1   1.00
      0.77904756001          1.0000000        
S   1   1.00
      0.30229502043          1.0000000        
P   1   1.00
     99.782996032            1.0000000
P   1   1.00                          
     23.176124101            1.0000000
P   1   1.00                          
      7.1163945872           1.0000000
P   1   1.00                          
      2.4418711435           1.0000000
P   1   1.00                          
      0.83389605766          1.0000000        
P   1   1.00
      0.26607311301          1.0000000        
D   1   1.00
      4.01400000             1.0000000        
D   1   1.00
      1.09600000             1.0000000        
F   1   1.00
      2.54400000             1.0000000        
****
H     0 
S   3   1.00
     34.0613410              0.60251978E-02   
      5.1235746              0.45021094E-01   
      1.1646626              0.20189726       
S   1   1.00
      0.32723041             1.0000000        
S   1   1.00
      0.10307241             1.0000000        
P   1   1.00
      0.8000000              1.0000000        
****
EOF

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
