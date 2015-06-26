#!/usr/bin/perl

use strict;
use warnings;
use POSIX qw/ceil/;
use POSIX qw/floor/;

open(OUT,"<",$ARGV[0]) or die $!;
my $FLAG=0;
my $Et;
my $OCC;
my $UNC;
my $OCCline;
my $LUMOFlag=0;
my $HOMO;
my $LUMO;

while(my $line=<OUT>){
  if($line =~ m/NUMBER OF ELECTRONS/){
    $Et = $line;
    $Et =~ s/.*: *([0-9]*)\..*/$1/g;
    chop $Et;
    $OCC = ceil($Et/2);
    $UNC = $OCC + 1;
    $OCCline = ceil($OCC/2);
    if(($OCC % 2)==0){
      $LUMOFlag=1;
    }
  }
  elsif($line =~ m/EIGENVALUES\(EV\) AND OCCUPATION/){
    for(my $i=0;$i<=$OCCline;$i++){
      $line = <OUT>;
    }
    $HOMO = $line;
    $HOMO =~ s/.*$OCC *([0-9\.]*) *.*/$1/;
    chop $HOMO;
    if($LUMOFlag==1){$line=<OUT>}
    $LUMO = $line;
    $LUMO =~ s/.*$UNC *([0-9\.]*) *.*/$1/;
    chop $LUMO;
  }
}
my $Eg = $LUMO - $HOMO;
print "$LUMO $HOMO\n";
