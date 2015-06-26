#!/usr/bin/perl
#
# USAGE: ./alignX.pl *.xyz A B C
#  aligning BA with positive x axis
#  while BC lies in xy plane
#


# check .xyz file exist
unless(-e $ARGV[0]){
  print "xyz file '$ARGV[0]' not found.\n\n";
  print "USAGE: ./alignX.pl toConvert.xyz\n";
  exit;
}

# open and read aligning coordinates
# read in number of atom $NA and 
# second line information $info2
open XYZ, "<", "$ARGV[0]" or die $!;
my $NA=<XYZ>;
$NA =~ s/[ \t\n]*//g;
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

# print result as xyz format
print "$NA\n";
print $info2;
for(my $a=0;$a<scalar(@ATOMS);$a++){
  printf("%-2s", $ATOMS[$a][0]);
  for(my $i=1;$i<4;$i++){printf(" % 11.8f",$ATOMS[$a][$i]*0.52917721092);}
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
