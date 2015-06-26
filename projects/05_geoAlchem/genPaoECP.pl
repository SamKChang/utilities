#!/usr/bin/perl

use warnings;
use strict;
use POSIX;

if($#ARGV != 1){
  print "two Gaussian output files expected\n";
  exit;
}

#if($#ARGV == 0){ 
#  unless(($ARGV[0] =~ /\.log/)||($ARGV[0] =~ /\.out/)){
#    print "foo.out or bar.log is expected.\n";
#    exit;
#  }
#}else{
#  print "foo.out and bar.log are expected\n"; 
#  exit;
#}
open my $gout, "<", $ARGV[0] or die $!;
open my $rout, "<", $ARGV[1] or die $!;

my $tarName=$ARGV[0];
$tarName =~ s/\..*//;

my ($NBasisT, $NBasisR,$NPart);
$NBasisT = $NPart = 0;
my (@ECPAOt, @ECPAOr);
my (@Vt,@Vr,$Prt,$Pt,$Pr);

while(my $line = <$gout>){
  if($line =~ m/Calculate potential energy integrals/){
    $NBasisT = <$gout>;
    $NBasisT =~ s/.*NBasis = *//;
    $NBasisT =~ s/\ .*\n//;
    $NPart = ceil($NBasisT/5);
    print "target basis: $NBasisT\n";
    for(my $i=0;$i<$NBasisT**2;$i++){
      $ECPAOr[$i]=0;
    }
    print "$NBasisT x $NBasisT target AO array initialized\n";
  }
  elsif($line =~ m/Potential Energy/){
    for(my $n=0;$n<$NPart;$n++){
      $line = <$gout>;
      for(my $i=$n*5;$i<$NBasisT;$i++){
        $line = <$gout>;
        $line =~ s/D/E/g;
        my @AOline = split(" ", $line);
        for(my $j=1;$j<6;$j++){
          if($j<scalar(@AOline)){
            my $s=$i + $NBasisT*($j-1 + $n*5);
            my $sj = int($s/$NBasisT);
            my $si = $s - $NBasisT*$sj;
            my $t = $sj + $NBasisT*$si;
            if($s==$t){
              $Vt[$s] = $AOline[$j];
            }
            else{
              $Vt[$s] = $AOline[$j];
              $Vt[$t] = $AOline[$j];
            }
          }
        }
      }
    }
  }
  elsif($line =~ m/ECP Integrals/){
    for(my $n=0;$n<$NPart;$n++){
      $line = <$gout>;
      for(my $i=$n*5;$i<$NBasisT;$i++){
        $line = <$gout>;
        $line =~ s/D/E/g;
        my @AOline = split(" ", $line);
        for(my $j=1;$j<6;$j++){
          if($j<scalar(@AOline)){
            my $s=$i + $NBasisT*($j-1 + $n*5);
            my $sj = int($s/$NBasisT);
            my $si = $s - $NBasisT*$sj;
            my $t = $sj + $NBasisT*$si;
            if($s==$t){
              $ECPAOt[$s] = $AOline[$j];
            }   
            else{
              $ECPAOt[$s] = $AOline[$j];
              $ECPAOt[$t] = $AOline[$j];
            }   
          }   
        }   
      }   
    }   
  }
}

print "\n";

while(my $line = <$rout>){
  if($line =~ m/Calculate potential energy integrals/){
    $NBasisR = <$rout>;
    $NBasisR =~ s/.*NBasis = *//;
    $NBasisR =~ s/\ .*\n//;
    $NPart = ceil($NBasisR/5);
    print "reference basis: $NBasisR\n";
    for(my $i=0;$i<$NBasisR**2;$i++){
      $ECPAOr[$i]=0;
    }
    print "$NBasisR x $NBasisR reference AO array initialized\n";
  }
  elsif($line =~ m/Potential Energy/){
    for(my $n=0;$n<$NPart;$n++){
      $line = <$rout>;
      for(my $i=$n*5;$i<$NBasisR;$i++){
        $line = <$rout>;
        $line =~ s/D/E/g;
        my @AOline = split(" ", $line);
        for(my $j=1;$j<6;$j++){
          if($j<scalar(@AOline)){
            my $s=$i + $NBasisR*($j-1 + $n*5);
            my $sj = int($s/$NBasisR);
            my $si = $s - $NBasisR*$sj;
            my $t = $sj + $NBasisR*$si;
            if($s==$t){
              $Vr[$s] = $AOline[$j];
            }
            else{
              $Vr[$s] = $AOline[$j];
              $Vr[$t] = $AOline[$j];
            }
          }
        }
      }
    }
  }
  elsif($line =~ m/ECP Integrals/){
    for(my $n=0;$n<$NPart;$n++){
      $line = <$rout>;
      for(my $i=$n*5;$i<$NBasisR;$i++){
        $line = <$rout>;
        $line =~ s/D/E/g;
        my @AOline = split(" ", $line);
        for(my $j=1;$j<6;$j++){
          if($j<scalar(@AOline)){
            my $s=$i + $NBasisR*($j-1 + $n*5);
            my $sj = int($s/$NBasisR);
            my $si = $s - $NBasisR*$sj;
            my $t = $sj + $NBasisR*$si;
            if($s==$t){
              $ECPAOr[$s] = $AOline[$j];
            }
            else{
              $ECPAOr[$s] = $AOline[$j];
              $ECPAOr[$t] = $AOline[$j];
            }
          }
        }
      }
    }
  }
}
open(my $fao, "> ./$tarName.pao") || die "open file\n";
for(my $j=0;$j<$NBasisR;$j++){
  for(my $i=0;$i<$NBasisR;$i++){
    my $sr = $i + $NBasisR*$j;
    my $st = $i + $NBasisT*$j;
    $Prt = -$Vt[$st] + $ECPAOt[$st];
    $Pr = -$Vr[$sr] + $ECPAOr[$sr];
    # perturbation: (Ht - Hr)
    $Pt =  $Prt - 2*$Pr;
    printf($fao "% 12.6E\n", $Pt);
  }
}

print "\n$NBasisR x $NBasisR (Vt - Vr) perturbation matrix generated\n";

print "output file: $tarName.pao\n";
