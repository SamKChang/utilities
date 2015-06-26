#!/usr/bin/perl

use strict;
use warnings;
use Chemistry::Mol;
use Chemistry::File::XYZ;
use Chemistry::Bond::Find ':all';

my $mol = Chemistry::Mol->read("$ARGV[0]");
my $stem = $ARGV[0];
$stem =~ s/\.xyz//;

find_bonds($mol);
assign_bond_orders($mol);

my @MOLS = $mol->separate;
for(my $i=0;$i<scalar(@MOLS);$i++){
  my $moli = $MOLS[$i];
  my $a = sprintf("%02d",$i+1);
  $moli->write("$stem-$a.xyz");
}
