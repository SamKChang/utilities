#!/usr/bin/perl
#
# USAGE: ./stretch.pl *.xyz A B length C D ....
#  stretching axix    : AB
#  centering atom     : A
#  stretching distance: length
#  stretching atoms   : C D ...

# check .xyz file exist
unless(-e $ARGV[0]){
  print "xyz file '$ARGV[0]' not found.\n\n";
  print "USAGE: ./stretch.pl *.xyz A B length C D ...\n";
  print "stretching axix    : AB\n";
  print "centering atom     : A\n";
  print "stretching distance: length\n";
  print "stretching atoms   : C D ...\n";
  exit;
}

# check input format
# number of moving atoms
my $N = @ARGV - 4;
for(my $i=1;$i<3;$i++){
  unless($ARGV[$i] =~ /^[1-9]\d*$/ ){
    print "$i"."^th atom index '$ARGV[$i]' is not a positive integer.\n";
    exit;
  }
}
for(my $i=4;$i<@ARGV;$i++){
  unless($ARGV[$i] =~ /^[1-9]\d*$/ ){
    my $j = $i-4;
    print "$j"."^th atom index '$ARGV[$i]' is not a positive integer.\n";
    exit;
  }
}

# read in number of atom $NA and 
# second line information $info2
open XYZ, "<", "$ARGV[0]" or die $!;
my $NA = <XYZ>;
my $info2 = <XYZ>;

# read atom type and coordinates
for(my $itr=0;$itr<$NA;$itr++){
  if($itr>=0){
    # read and string into componets
    my @ATOM = split(' ',<XYZ>);
    # assign multiD data array
    push (@{$ATOMS[$itr]}, @ATOM);
  }
}

my @infoEND;
for($itr=0;<XYZ>;$itr++){
  #my @info = split(' ',$_);
  push (@infoEND, $_);
}
my $Ni=$itr;

my @vecAB;
my $lAB = 0;
for(my $i=0;$i<3;$i++){
  $vecAB[$i] = $ATOMS[$ARGV[2]-1][$i+1] - $ATOMS[$ARGV[1]-1][$i+1];
  $lAB += $vecAB[$i]**2;
}
$lAB = $ARGV[3] * $lAB**(-1/2);

# rotate atoms w.r AB axis
my @Ncoord = map {[@$_]} @ATOMS;
for(my $a=0;$a<$N;$a++){
  my $aa = $ARGV[4+$a] - 1;
  for(my $i=0;$i<3;$i++){
    #$Ncoord[$aa][$i+1] = $ATOMS[$aa][$i+1] + $vecAB[$i]*$lAB;
    $Ncoord[$aa][$i+1] = $ATOMS[$aa][$i+1] -$ATOMS[$ARGV[2]-1][$i+1] + $vecAB[$i]*$lAB;
  }
}

# print result as xyz format
print $NA;
print $info2;
for(my $a=0;$a<scalar(@ATOMS);$a++){
  printf("%-2s", $Ncoord[$a][0]);
  for(my $i=1;$i<4;$i++){printf(" % 11.8f",$Ncoord[$a][$i]);}
  print "\n";
}
for($itr=0;$itr<$Ni;$itr++){
  print "$infoEND[$itr]";
}

close(XYZ);
