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

# read atom type and coordinates
for(my $itr=0;$itr<$NA;$itr++){
  # read and string into componets
  my @ATOM = split(' ',<XYZ>);
  # assign multiD data array
  push (@{$ATOMS[$itr]}, @ATOM);
  $Nve += atom_ve($ATOM[0]);
}
$multi += $Nve%2;

# set multiplicity to 1 by adding electrons
if(($ARGV[0] =~ m/^cat/)&&($multi != 1)){
  $multi = 1;
  $charge = 1;
}elsif($multi != 1){
  $multi = 1;
  $charge = -1;
}

# read what every info at the end of file
my @infoEND;
for($itr=0;<XYZ>;$itr++){
  #my @info = split(' ',$_);
  push (@infoEND, $_);
}
my $Ni=$itr;

my $pwd = cwd();



#############################
# print nwchem input format #
#############################
print <<EOF;
start $stem
title $stem

geometry units an nocenter noautoz
EOF


# print atom type and coordinates #
for(my $a=0;$a<$NA;$a++){
  printf(" %-2s", $ATOMS[$a][0]);
  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS[$a][$i]);}
  print "\n";
}
print "end\n";

if($charge != 0){
  print "charge $charge\n";
}


# print basis format
print "basis cartesian\n";
#for(my $a=0;$a<$NA;$a++){
#  print(" %-2s library Def2-QZVPD\n");
#}


# print directive
print <<EOF;
end
dft
 xc xpbe96 cpbe96
 vectors output $stem.movecs
end
task dft
EOF

# density output
#dplot
# vectors scr$stem.movecs
#  LimitXYZ
#  -5 8 150
#  -5 5 100
#  -5 5 100
# gaussian
# spin total
# output $stem.cube
#end
#task dplot

#print "\n";
## extract 0th column from ATOMS
## map $_->[0], @ATOMS
#my @col = uniq(map $_->[0], @ATOMS);
#foreach(@col){
#  if($_ =~ m/H/){
#    atom_HBasis($_);
#  }
#  else{
#    atom_NeBasisU($_);
#  }
#}

close(XYZ);

###################
# useful routines #
###################
# uniq function
sub uniq {
  my %temp_hash = map { $_ => 1 } @_; 
  return keys %temp_hash;
}

# count valence electrons
sub atom_ve{
  switch($_[0]){
    case "H" {return 1}
    case "Li" {return 1}
    case "Be" {return 2}
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

sub atom_a2n{
  switch($_[0]){
    case "H" {return  1}
    case "Li"{return  3}
    case "Be"{return  4}
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
    case "3" {return "Li" }
    case "4" {return "Be" }
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
    else    {
      print "element: $_[0] not found\n";
      die;
    }
  }
}
