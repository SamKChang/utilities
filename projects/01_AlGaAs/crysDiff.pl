#!/usr/bin/perl

use strict;
use warnings;

sub readCoord($);
my $cut=10**-8;

my ($coord, $na, $elem, $pps)=readCoord($ARGV[0]);
my @coord1=@$coord;
my @na1=@$na;
my @elem1=@$elem;
my @PPs1=@$pps;
($coord, $na, $elem, $pps)=readCoord($ARGV[1]);
my @coord2=@$coord;
my @na2=@$na;
my @elem2=@$elem;

my @coord3=@coord1;
my @na3=@na1;
my @elem3=@elem1;
my @list;

# loop through all possible atom pair
for(my $N1=0;$N1<scalar @na1;$N1++){
  for(my $N2=0;$N2<scalar @na2;$N2++){
    for(my $I=0;$I<scalar $na1[$N1];$I++){
      for(my $J=0;$J<scalar $na2[$N2];$J++){
        # pair wise distance between two atoms
        my $dR=0;
        for(my $i=0;$i<3;$i++){
          $dR+=abs($coord1[$N1][$I][$i]-$coord2[$N2][$J][$i]);
        }
        # if coordinate matched, proceed with pairing
        if($dR<=$cut and $N1!=$N2){
          my @coord=@{$coord1[$N1][$I]};
          my $type="$elem1[$N1]2$elem2[$N2]";
          # check whether $type is already in element list
          my @index=grep {$elem3[$_] =~ /$type/} 0..$#elem3;
          if(scalar @index==0){
            my $N3=scalar @coord3;
            push @elem3, $type;
            push @na3, 1;
            push @{$coord3[$N3][0]}, @coord;
          }else{
            my $I3=$index[0];
            my $i3=scalar @{$coord3[$I3]};
            $na3[$I3]++;
            push @{$coord3[$I3][$i3]}, @coord;
          }
          # list of coordinates to remove
          push @list, [$N1,$I];
          $na3[$N1]--;

          # print matched coordinates
          #print "$N1 $N2, $I $J: ";
          #print "$elem1[$N1]2$elem2[$N2]->";
          #for(my $i=0;$i<3;$i++){
          #  print("  $coord1[$N1][$I][$i]");
          #}
          #print "\n";
        }
      }
    }
  }
}

# sort multidimensional array by desending order
@list = sort {$b->[1] <=> $a->[1]} @list;
for(my $i=0;$i<scalar @list;$i++){
  # remove element from array by index
  splice @{$coord3[$list[$i][0]]}, $list[$i][1], 1;
}

# print cpmd input format
open FILE, "<", $ARGV[0];
while(<FILE>){
  unless($_ =~ m/&ATOMS/){
    print $_;
  }else{
    last;
  }
}
close FILE;

# print atom type and coordinates
print "&ATOMS\n";
for(my $N3=0;$N3<scalar @na3;$N3++){
  if($N3<@na1){
    print "$PPs1[$N3] LMAX=F\n  $na3[$N3]\n";
  }else{
    print "*$elem3[$N3]_000.psp FRAC\n LMAX=F\n  $na3[$N3]\n";
  }
  for(my $I=0;$I<$na3[$N3];$I++){
    print " ";
    for(my $i=0;$i<3;$i++){
      print "  $coord3[$N3][$I][$i]"
    }
    print "\n";
  }
  print "\n";
}
print "&END\n";

sub readCoord($) {
  open FILE, "<", $_[0];
  my (@atoms, @Na, @COORDS, @Natom, @PPs);
  my $I=0;
  while(<FILE>){
    if($_ =~ m/^\*/){
      my $line = $_;
      push @PPs, $line;
      $line =~ s/^\*([A-Za-z]*)_.*/$1/;
      chop $line;
      push @atoms, $line;
      my $Na=<FILE>;
      $Na=<FILE>;
      $Na =~ s/ *//;
      for(my $i=0;$i<$Na;$i++){
        $line = <FILE>;
        my @coord = split " ", $line;
        push(@{$COORDS[$I][$i]},@coord);
      }
      push @Natom, $Na;
      $I++
    }
  }
  close FILE;
  return (\@COORDS,\@Natom,\@atoms, \@PPs);
}
