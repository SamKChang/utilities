#!/usr/bin/perl -w

#use strict;
use warnings;
use Switch;
use Math::SpecFun::Erf qw(erf);

# expect two input foo.out boo.fchk
# gaussian09 format are expected

if($#ARGV == 0){
  if($ARGV[0] !~ /\.fchk/){
    print "foo.out and bar.fchk are expected.\n";
    exit;
  }
}else{
  print "foo.out and bar.fchk are expected\n"; 
  exit;
}
#open my $gout, "<", $ARGV[0] or die $!;
open my $fchk, "<", $ARGV[0] or die $!;
$ARGV[0] =~ s/fchk/xyz/;

print "$ARGV[0]\n";

my $TAO;
my $multy=0;
my @CMNBasis;

my $acharge;
my @ACHARGE;
my $type;
my @TYPE;
my $nshell;
my @NSHELL;
my $map;
my @MAP;
my $pexp;
my @PEXP;
my $ctrc;
my @CTRC;
my $coord;
my @COORD;
my $ew_a;
my @EW_a;
my $ev_a;
my @EV_a;
my $ew_b;
my @EW_b;
my $ev_b;
my @EV_b;

my $NAO;
my $NPG;
my $NCRD;

my $NMO;
my $NMOl;
my $NATOM;

# read in from Gaussian09 foo.fchk file
my $Diff = 0;
my $preLine = <$fchk>;
while (my $line = <$fchk>){
  if($line =~ m/^Charge/){
    $line =~ s/.*I *//;
    $CMNBasis[0] = $line;
    for(my $i=1;$i<6;$i++){
      $line = <$fchk>;
      $line =~ s/.*I *//;
      $CMNBasis[$i] = $line;
    }
    $multy = $CMNBasis[1];
  }elsif($preLine =~ m/^Nuclear charges/){
    $preLine =~ s/.*N= *//;
    unless($preLine % 5){$Diff = 1;}
    else{$Diff = 0;}
    my $N_tmp = int($preLine / 5) - $Diff;
    $acharge = $line;
    for(my $i=0;$i<$N_tmp;$i++){
      $acharge .= <$fchk>;
    }
    @ACHARGE = split(" ", $acharge);
  }elsif($preLine =~ m/^Shell types/){
    $preLine =~ s/.*N= *//;
    $TAO = $preLine;
    unless($TAO % 6){$Diff = 1;}
    else{$Diff = 0;}
    $NAO = int($preLine / 6) - $Diff;
    $type = $line;
    for(my $i=0;$i<$NAO;$i++){
      $type .= <$fchk>;
    }
    @TYPE = split(" ", $type)
  }elsif($preLine =~ m/^Number of primitives per shell/){
    $nshell = $line;
    for(my $i=0;$i<$NAO;$i++){
      $nshell .= <$fchk>;
    }
    @NSHELL = split(" ", $nshell)
  }elsif($preLine =~ m/^Shell to atom map/){
    $map = $line;
    for(my $i=0;$i<$NAO;$i++){
      $map .= <$fchk>;
    }
    @MAP = split(" ", $map)
  }elsif($preLine =~ m/^Primitive exponents/){
    $preLine =~ s/.*N= *//;
    unless($preLine % 5){$Diff = 1;}
    else{$Diff = 0;}
    $NPG = int($preLine / 5) - $Diff;
    $pexp = $line;
    for(my $i=0;$i<$NPG;$i++){
      $pexp .= <$fchk>;
    }
    @PEXP = split(" ", $pexp)
  }elsif($preLine =~ m/^Contraction coefficients/){
    $ctrc = $line;
    for(my $i=0;$i<$NPG;$i++){
      $ctrc .= <$fchk>;
    }
    @CTRC = split(" ", $ctrc);
  }elsif($preLine =~ m/^Current cartesian coordinates/){
    $preLine =~ s/.*N= *//;
    $NATOM = $preLine/3;
    unless($preLine % 5){$Diff = 1;}
    else{$Diff = 0;}
    $NCRD = int($preLine / 5) - $Diff;
    $coord = $line;
    for(my $i=0;$i<$NCRD;$i++){
      $coord .= <$fchk>;
    }
    @COORD = split(" ", $coord);
  }elsif($preLine =~ m/^Alpha Orbital Energies/){
    $preLine =~ s/.*N= *//;
    $NMO = $preLine;
    unless($NMO % 5){$Diff = 1;}
    else{$Diff = 0;}
    $NMOl= int($NMO / 5) - $Diff;
    $ew_a = $line;
    for(my $i=0;$i<$NMOl;$i++){
      $ew_a .= <$fchk>;
    }
    @EW_a = split(" ", $ew_a);
  }elsif($preLine =~ m/^Alpha MO coefficients/){
    unless($NMO**2 % 5){$Diff = 1;}
    else{$Diff = 0;}
    my $NMOAO= int($NMO**2 / 5) - $Diff;
    $ev_a = $line;
    for(my $i=0;$i<$NMOAO;$i++){
      $ev_a .= <$fchk>;
    }
    @EV_a = split(" ", $ev_a);
  }elsif(($preLine =~ m/^Beta MO coefficients/)&($multy!=1)){
    unless($NMO**2 % 5){$Diff = 1;}
    else{$Diff = 0;}
    my $NMOAO= int($NMO**2 / 5) - $Diff;
    $ev_b = $line;
    for(my $i=0;$i<$NMOAO;$i++){
      $ev_b .= <$fchk>;
    }
    @EV_b = split(" ", $ev_b);
  }
  $preLine = $line;
}

# list of AOs string
# not necessary
#my @ORBITALS;
#my @orbital;
#while(my $line = <$gout>){
#  if($line =~ m/ Molecular Orbital Coefficients:/){
#    for(my $i=0;$i<3;$i++){
#      $line = <$gout>;
#    }
#    $line = <$gout>;
#    @ORBITALS = $line =~ m/[0-9]{1,2}[SPXYZ]{1,4}/g;
#    for(my $i=1;$i<$NMO;$i++){
#      $line = <$gout>;
#      @orbital = $line =~ m/[0-9]{1,2}[SPXYZ]{1,4}/g;
#      push(@ORBITALS, "$orbital[0]");
#    }
#  }
#}

# sum up number of AOs
my $sum = 0;
my $j=1;
my @AOMAP;
my @lmMAP;
foreach(@NSHELL){
  switch($TYPE[$j-1]){
    case 0 {
      $sum += 1;
      push(@AOMAP, $j);
      push(@lmMAP, '0 0 0');
    }
    case 1 {
      $sum += 3;
      push(@AOMAP, $j,$j,$j);
      push(@lmMAP, '1 0 0', '0 1 0', '0 0 1');
    }
    case 2 {
      $sum += 6;
      push(@AOMAP, $j,$j,$j,$j,$j,$j);
      push(@lmMAP,'2 0 0','0 2 0','0 0 2');
      push(@lmMAP,'1 1 0','1 0 1','0 1 1');
    }
    case 3 {
      $sum += 10;
      push(@AOMAP, $j,$j,$j,$j,$j,$j,$j,$j,$j,$j);
      push(@lmMAP,'3 0 0','0 3 0','0 0 3','1 2 0','2 1 0');
      push(@lmMAP,'2 0 1','1 0 2','0 1 2','0 2 1','1 1 1');
    }
  }
  $j++;
}

# AGMAP: AO to primitive Gaussain
my $sum_i=0;
my @AGMAP=0;
foreach(@NSHELL){
  $sum_i += $_;
  push(@AGMAP,$sum_i);
}

#foreach(@AOMAP){
#  print $_,"\n";
#}

#for(my $i=0;$i<$NMO;$i++){
#  for(my $j=$i;$j<$NMO;$j++){
#    #if($lmMAP[$i] =~ $lmMAP[$j]){
#      print "<",$i+1,"|",$j+1,">:\t",$lmMAP[$i],"-",$lmMAP[$j],"\t";
#      my $AO_i=$AOMAP[$i]-1;
#      my $AO_j=$AOMAP[$j]-1;
#      if ($AGMAP[$AO_i]==$AGMAP[$AO_j]){
#        for(my $m=$AGMAP[$AO_i];$m<$AGMAP[$AO_i+1];$m++){
#          for(my $n=$m;$n<$AGMAP[$AO_j+1];$n++){
#            print " ($m|$n)";
#          }
#        }
#        print "\n";
#      }else{
#        for(my $m=$AGMAP[$AO_i];$m<$AGMAP[$AO_i+1];$m++){
#          for(my $n=$AGMAP[$AO_j];$n<$AGMAP[$AO_j+1];$n++){
#            print " ($m|$n)";
#          }
#        }
#        print "\n";
#      }
#    #}
#  }
#}

#my $c=0;
#open(my $fcmn, "> ./CMNBasis.pltmp") || die "open file\n";
#open(my $facharge, "> ./ACHARGE.pltmp") || die "open file\n";
#open(my $ftype, "> ./TYPE.pltmp") || die "open file\n";
#open(my $fnshell, "> ./NSHELL.pltmp") || die "open file\n";
#open(my $fmap, "> ./MAP.pltmp") || die "open file\n";
open(my $fcoord, "> ./$ARGV[0]") || die "open file\n";
#open(my $fpexp, "> ./PEXP.pltmp") || die "open file\n";
#open(my $fctrc, "> ./CTRC.pltmp") || die "open file\n";
#open(my $fewa, "> ./EW_a.pltmp") || die "open file\n";
#open(my $feva, "> ./EV_a.pltmp") || die "open file\n";
#open(my $faomap, "> ./AOMAP.pltmp") || die "open file\n";
#open(my $flmmap, "> ./lmMAP.pltmp") || die "open file\n";
#open(my $fagmap, "> ./AGMAP.pltmp") || die "open file\n";
#foreach(@CMNBasis){
#  print $fcmn "$_";
#}
#for(my $i=0;$i<scalar(@TYPE);$i++){
#  print $ftype "$TYPE[$i]\n";
#  print $fnshell "$NSHELL[$i]\n";
#  print $fmap "$MAP[$i]\n";
#  print $fagmap "$AGMAP[$i]\n";
#  $c += 3;
#}
#print $fagmap "$AGMAP[scalar(@TYPE)]";
#for(my $i=0;$i<scalar(@PEXP);$i++){
#  print $fpexp "$PEXP[$i]\n";
#  print $fctrc "$CTRC[$i]\n";
#}
#for(my $i=0;$i<scalar(@EW_a);$i++){
#  print $fewa "$EW_a[$i]\n";
#  print $faomap "$AOMAP[$i]\n";
#  print $flmmap "$lmMAP[$i]\n";
#}
#for(my $i=0;$i<scalar(@EV_a);$i++){
#  print $feva "$EV_a[$i]\n";
#}


my @ATYPE;
my $itype = 0;
foreach(@ACHARGE){
  switch(int($_)){
    case 1 {
      $ATYPE[$itype] = 'H';
    }
    case 6 {
      $ATYPE[$itype] = 'C';
    }
    case 7 {
      $ATYPE[$itype] = 'N';
    }
    case 8 {
      $ATYPE[$itype] = 'O';
   }
  }
  $itype++;
}
my $a = 0.52917721092;
printf $fcoord "%d\n\n", $NATOM;
for(my $i=0;$i<$NATOM;$i++){
  printf $fcoord "%s % 14.8E % 14.8E % 14.8E\n", $ATYPE[$i], $COORD[3*$i]*$a, $COORD[3*$i+1]*$a, $COORD[3*$i+2]*$a;
}
#  print $fcoord "$ACHARGE[$i]\n";
#  print $fcoord "$COORD[3*$i] ";
#  print $fcoord "$COORD[3*$i+1] ";
#  print $fcoord "$COORD[3*$i+2]\n";
#if($multy!=1){
#  open(my $fewb, "> ./EW_b.pltmp") || die "open file\n";
#  open(my $fevb, "> ./EV_b.pltmp") || die "open file\n";
#  foreach(@EW_b){
#    print $fewb "$_\n";
#  }
#  foreach(@EV_b){
#    print $fevb "$_\n";
#  }
#}

#
#my $J=0;
#foreach(@AOMAP){
#  print "$lmMAP[$J]\t$_\n";
#  $J++;
#}

#my $j=0;
#my $ao=0;
#foreach(@NSHELL){
#  switch($TYPE[$j]){
#    case 0 {
#      print ++$ao, ": S\n";
#    }
#    case 1 {
#      print ++$ao, ": P\n";
#      print ++$ao, ": P\n";
#      print ++$ao, ": P\n";
#    }                  
#    case 2 {           
#      print ++$ao, ": D\n";
#      print ++$ao, ": D\n";
#      print ++$ao, ": D\n";
#      print ++$ao, ": D\n";
#      print ++$ao, ": D\n";
#      print ++$ao, ": D\n";
#    }                  
#    case 3 {           
#      print ++$ao, ": F\n";
#      print ++$ao, ": F\n";
#      print ++$ao, ": F\n";
#      print ++$ao, ": F\n";
#      print ++$ao, ": F\n";
#      print ++$ao, ": F\n";
#      print ++$ao, ": F\n";
#      print ++$ao, ": F\n";
#      print ++$ao, ": F\n";
#      print ++$ao, ": F\n";
#    }
#  }
#  $j++;
#}

#my $j=0;
#my $g=0;
#foreach (@NSHELL){
#  my $x = $j*3;
#  print "$TYPE[$j] $NSHELL[$j] $MAP[$j]: (";
#  for(my $i=0;$i<3;$i++){
#    print " $COORD[$x+$i]";
#  }
#  print "), gaussians: (";
#  for(my $s=0;$s<$NSHELL[$j];$s++){
#    print " $CTRC[$g]*exp($PEXP[$g])";
#    $g = $g + 1;
#  }
#  print ")\n";
#  $j++;
#}
