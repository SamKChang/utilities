#!/usr/bin/perl

use strict;
use warnings;
use Chemistry::Mol;
use Chemistry::File::XYZ;
use Chemistry::Bond::Find ':all';
use Data::Dumper;

my $mol1 = Chemistry::Mol->read("$ARGV[0]");
my $mol2 = Chemistry::Mol->read("$ARGV[1]");
my $stem1 = $ARGV[0];
$stem1 =~ s/\.xyz//;
my $stem2 = $ARGV[1];
$stem2 =~ s/\.xyz//;

find_bonds($mol1);
assign_bond_orders($mol1);
find_bonds($mol2);
assign_bond_orders($mol2);

my @HBondDA1 = grep {
  $_->symbol =~ m/[ONS]/ 
} $mol1->atoms;

my @HBondDA2 = grep {
  $_->symbol =~ m/[ONS]/ 
} $mol2->atoms;

# number of H-bonds
my $NHBonds = 0;

# H-bond doners from mol1
for(my $i=0;$i<scalar(@HBondDA1);$i++){
  # H-bond acceptors from mol2
  for(my $j=0;$j<scalar(@HBondDA2);$j++){
    # check doner-acceptor distance
    if($HBondDA1[$i]->distance($HBondDA2[$j])<3.5){
      # check doner-hydrogen-acceptor distance
      my @N1 = grep {
        $_->symbol =~ m/H/
        and $_->distance($HBondDA2[$j]) < 3.0
      } $HBondDA1[$i]->neighbors;
      my @N2 = grep {
        $_->symbol =~ m/H/
        and $_->distance($HBondDA1[$i]) < 3.0
      } $HBondDA2[$j]->neighbors;
      $NHBonds += (scalar(@N1)+scalar(@N2));
    }
  }
}

if($NHBonds > 0){
  my $mol3 = $mol1->combine($mol2);
  print "$stem1-$stem2: $NHBonds\n";
  $mol3->write("$stem1-H-$stem2.xyz");
}



#foreach(@close_hydrogens){
#  print $_->symbol," ", $_->coords,"\n";
#}
#
#foreach($mol->bonds){
#  print $_,"\n";
#}
#
#print $mol->atoms(1)->distance($mol->atoms(4));
#
#no strict 'refs';
#print "Instance METHOD IS  " . Dumper( \%{ref ($mol->atoms(1))."::" });
#
#print $mol->atoms(150)->angle_deg(148,82), "\n";
#
#
#print "\n";
#
##my @MOLS = $mol->separate;
##for(my $i=0;$i<scalar(@MOLS);$i++){
##  my $moli = $MOLS[$i];
##  my $a = sprintf("%02d",$i+1);
##  $moli->write("$stem-$a.xyz");
##}
