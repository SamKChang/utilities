#!/usr/bin/perl
#
# USAGE: ./alignX.pl *.xyz A B
#  aligning BA with positive x axis
#  while BC lies in xy plane
#


# check .xyz file exist
unless(-e $ARGV[0]){
  print "xyz file '$ARGV[0]' not found.\n\n";
  print "USAGE: ./alignXlin.pl *.xyz A B\n";
  print " aligning AB with positive x axis\n";
  exit;
}

# check posititive number of atom list
for(my $i=1;$i<2;$i++){
  unless($ARGV[$i] =~ /^[1-9]\d*$/ ){
    print "$i"."^th atom index '$ARGV[$i]' is not a positive integer.\n";
    exit;
  }
}


# open and read aligning coordinates
# read in number of atom $NA and 
# second line information $info2
open XYZ, "<", "$ARGV[0]" or die $!;
my $NA=<XYZ>;
my $info2=<XYZ>;
for(my $itr=0;$itr<$NA;$itr++){
  # read and string into componets
  my @ATOM = split(' ',<XYZ>);
  # assign multiD data array
  push (@{$ATOMS[$itr]}, @ATOM);
}
my @infoEND;
for($itr=0;<XYZ>;$itr++){
  #my @info = split(' ',$_);
  push (@infoEND, $_);
}
my $Ni=$itr;

# centering at A atom
my @Ccoord1 = map{[@$_]} @ATOMS;
for(my $a=0;$a<scalar(@ATOMS);$a++){
  for(my $i=1;$i<4;$i++){
    $Ccoord1[$a][$i] = $ATOMS[$a][$i] - $ATOMS[$ARGV[1]-1][$i];
  }
}

# rotate AB to XY plane
my $a1 = atan2($Ccoord1[$ARGV[2]-1][3],$Ccoord1[$ARGV[2]-1][1]);
my @R = GetRy($a1);
my @Rcoord1 = map {[@$_]} @Ccoord1;
for(my $a=0;$a<scalar(@ATOMS);$a++){
  for(my $i=0;$i<3;$i++){
    my $sum_i = 0;
    for(my $j=0;$j<3;$j++){
      $sum_i = $sum_i + $R[$i][$j] * $Ccoord1[$a][$j+1];
    }
    $Rcoord1[$a][$i+1] = $sum_i;
  }
}

# rotate AB to X axis
my $a1 = atan2($Rcoord1[$ARGV[2]-1][2],$Rcoord1[$ARGV[2]-1][1]);
my @R = GetRz(-$a1);
my @Rcoord2 = map {[@$_]} @Rcoord1;
for(my $a=0;$a<scalar(@ATOMS);$a++){
  for(my $i=0;$i<3;$i++){
    my $sum_i = 0;
    for(my $j=0;$j<3;$j++){
      $sum_i = $sum_i + $R[$i][$j] * $Rcoord1[$a][$j+1];
    }
    $Rcoord2[$a][$i+1] = $sum_i;
  }
}

# print result as xyz format
print $NA;
print $info2;
for(my $a=0;$a<scalar(@ATOMS);$a++){
  printf("%-2s", $Rcoord2[$a][0]);
  for(my $i=1;$i<4;$i++){printf(" % 11.8f",$Rcoord2[$a][$i]);}
  print "\n";
}
for($itr=0;$itr<$Ni;$itr++){
  print "$infoEND[$itr]";
}


## print result as xyz format
#print scalar(@ATOMS)."\n\n";
#for(my $a=0;$a<scalar(@ATOMS);$a++){
#  print $Rcoord3[$a][0];
#  #for(my $i=1;$i<4;$i++){printf("  % 14.8E",$Rcoord3[$a][$i]);}
#  for(my $i=1;$i<4;$i++){printf("  % 10.7f",$Rcoord3[$a][$i]);}
#  print "\n";
#}

#####################
# Rotation matrices #
#####################
#
sub GetRx {
  my $angle = $_[0];
  my $Rout = [  ];

  $Rout[0][0] = 1;
  $Rout[0][1] = 0;
  $Rout[0][2] = 0;

  $Rout[1][0] = 0;
  $Rout[1][1] = cos($angle);
  $Rout[1][2] = -sin($angle);

  $Rout[2][0] = 0;
  $Rout[2][1] = sin($angle);
  $Rout[2][2] = cos($angle);

  return @Rout;
}


sub GetRy {
  my $angle = $_[0];
  my $Rout = [  ];

  $Rout[0][0] = cos($angle);
  $Rout[0][1] = 0;
  $Rout[0][2] = sin($angle);
  
  $Rout[1][0] = 0;
  $Rout[1][1] = 1;
  $Rout[1][2] = 0;
  
  $Rout[2][0] = -sin($angle);
  $Rout[2][1] = 0;
  $Rout[2][2] = cos($angle);

  return @Rout;
}

sub GetRz {
  my $angle = $_[0];
  my $Rout = [  ];  

  $Rout[0][0] = cos($angle);
  $Rout[0][1] = -sin($angle);
  $Rout[0][2] = 0;
  
  $Rout[1][0] = sin($angle);
  $Rout[1][1] = cos($angle);
  $Rout[1][2] = 0;

  $Rout[2][0] = 0;
  $Rout[2][1] = 0;
  $Rout[2][2] = 1;
  
  return @Rout;
}

close(XYZ);