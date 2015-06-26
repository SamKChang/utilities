#!/usr/bin/perl
#
# USAGE: ./rotate.pl *.xyz A B deg C D ....
#  rotation axix:   AB
#  rotation center: A
#  rotating atoms:  C D ...

# check .xyz file exist
unless(-e $ARGV[0]){
  print "xyz file '$ARGV[0]' not found.\n\n";
  print "USAGE: ./rotate.pl *.xyz A B deg C D ....\n";
  print " rotation axix:   AB\n";
  print " rotation center: A\n";
  print " rotating atoms:  C D ...\n";
  exit;
}

# check input format
my $N = @ARGV - 4;
for(my $i=1;$i<3;$i++){
  unless($ARGV[$i] =~ /^[1-9]\d*$/ ){
    print "$i"."^th atom index '$ARGV[$i]' for rotation axix\n is not a positive integer.\n";
    exit;
  }
}

unless(($ARGV[3] <= 360)&&($ARGV[3] >= 0)){
  print "rotation angle: 0 < $ARGV[3] < 360 not valid\n";
  exit;
}


for(my $i=4;$i<@ARGV;$i++){
  unless($ARGV[$i] =~ /^[1-9]\d*$/ ){
    my $j = $i-4;
    print "$j"."^th index '$ARGV[$i]' for rotating atom is not a positive integer.\n";
    exit;
  }
}
 
# open and read aligning coordinates
# read in number of atom $NA and
# second line information $info2
open XYZ, "<", "$ARGV[0]" or die $!;
my $NA = <XYZ>;
my $info2 = <XYZ>;
for(my $itr=0;$itr<$NA;$itr++){
  if($itr>=0){
    # read and string into componets
    my @ATOM = split(' ',<XYZ>);
    # assign multiD data array
    push (@{$ATOMS[$itr]}, @ATOM);
  }
}

my $infoEND;
for($itr=0;<XYZ>;$itr++){
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
my $a1y = atan2($Ccoord1[$ARGV[2]-1][3],$Ccoord1[$ARGV[2]-1][1]);
my @R = GetRy($a1y);
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
my $a2z = atan2($Rcoord1[$ARGV[2]-1][2],$Rcoord1[$ARGV[2]-1][1]);
my @R = GetRz(-$a2z);
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
 
# rotate atoms w.r AB axis
my $a3x = 3.1415926535 * $ARGV[3] / 180;
my @R = GetRx($a3x);
my @Rcoord3 = map {[@$_]} @Rcoord2;
for(my $a=0;$a<$N;$a++){
  my $aa = $ARGV[4+$a] - 1;
  for(my $i=0;$i<3;$i++){
    my $sum_i = 0;
    for(my $j=0;$j<3;$j++){
      $sum_i = $sum_i + $R[$i][$j] * $Rcoord2[$aa][$j+1];
    }
    $Rcoord3[$aa][$i+1] = $sum_i;
  }
}

# rotate AB back to XY plane
my @R = GetRz($a2z);
my @Rcoord2 = map {[@$_]} @Rcoord3;
for(my $a=0;$a<scalar(@ATOMS);$a++){
  for(my $i=0;$i<3;$i++){
    my $sum_i = 0;
    for(my $j=0;$j<3;$j++){
      $sum_i = $sum_i + $R[$i][$j] * $Rcoord3[$a][$j+1];
    }   
    $Rcoord2[$a][$i+1] = $sum_i;
  }
}

# rotate AB back to original direction
my @R = GetRy(-$a1y);
my @Rcoord1 = map {[@$_]} @Rcoord2;
for(my $a=0;$a<scalar(@ATOMS);$a++){
  for(my $i=0;$i<3;$i++){
    my $sum_i = 0;
    for(my $j=0;$j<3;$j++){
      $sum_i = $sum_i + $R[$i][$j] * $Rcoord2[$a][$j+1];
    }   
    $Rcoord1[$a][$i+1] = $sum_i;
  }
}

# shift A back to old position
my @Ccoord1 = map{[@$_]} @Rcoord1;
for(my $a=0;$a<scalar(@ATOMS);$a++){
  for(my $i=1;$i<4;$i++){
    $Ccoord1[$a][$i] = $Rcoord1[$a][$i] + $ATOMS[$ARGV[1]-1][$i];
  }
}


# print result as xyz format
print scalar(@ATOMS)."\n";
print $info2;
for(my $a=0;$a<$NA;$a++){
  print $Rcoord3[$a][0];
  for(my $i=1;$i<4;$i++){printf("  % 11.8f",$Ccoord1[$a][$i]);}
  print "\n";
}
for($itr=0;$itr<$Ni;$itr++){
  print "$infoEND[$itr]";
}

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
