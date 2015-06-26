#!/usr/bin/perl

use strict;
use warnings;
use Chemistry::Mol;
use Chemistry::File::XYZ;
use Chemistry::Bond::Find ':all';
use Data::Dumper;

my $mol = Chemistry::Mol->read("$ARGV[0]");
my $stem = $ARGV[0];
$stem =~ s/\.xyz//;

find_bonds($mol);
#assign_bond_orders($mol);

print $mol;


my @close_hydrogens = grep {
  $_->symbol =~ m/[HC]/ 
  and $_->distance($mol->atoms(130)) < 3.0
} $mol->atoms;

foreach(@close_hydrogens){
  print $_->symbol," ", $_->coords,"\n";
}

foreach($mol->bonds){
  print $_,"\n";
}

print $mol->atoms(1)->distance($mol->atoms(4));

no strict 'refs';
print "Instance METHOD IS  " . Dumper( \%{ref ($mol)."::" });

print $mol->atoms(150)->angle_deg(148,82), "\n";


print "\n";

#my @MOLS = $mol->separate;
#for(my $i=0;$i<scalar(@MOLS);$i++){
#  my $moli = $MOLS[$i];
#  my $a = sprintf("%02d",$i+1);
#  $moli->write("$stem-$a.xyz");
#}
