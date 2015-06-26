#!/usr/bin/perl

use strict;
use warnings;

use WWW::Mechanize;

#my $url = "http://www.tcm.phy.cam.ac.uk/~mdt26/pseudo_lib/si/pseudo.html";
my $url = $ARGV[0];

my $mech = WWW::Mechanize->new();

$mech->get( $url );

my @links = $mech->links();

foreach my $link (@links) {
  my $line = $link->url();
  $line =~ s/.*\///;
  if($line =~ m/awfn/){
    print "$line\n";
  }
}
