#!/usr/bin/perl -w

#use strict;
use warnings;
use Switch;

# expect two input foo.out boo.fchk
# gaussian09 format are expected

if($#ARGV == 0){
  if($ARGV[0] !~ /\.fchk/){
    print "foo.fchk is expected.\n";
    print "usage: genCrd.pl foo.fchk\n";
    exit;
  }
}else{
  print "foo.fchk is expected\n"; 
  print "usage: genCrd.pl foo.fchk\n";
  exit;
}
#open my $gout, "<", $ARGV[0] or die $!;
open my $fchk, "<", $ARGV[0] or die $!;
my $out = $ARGV[0];
$out =~ s/\.fchk//;
my $crdOut = $out . ".crd";
my $pplOut = $out . ".ppl";
my $alkOut = $out . ".alk";
my $expOut = $out . ".exp";
my $kiOut  = $out . ".ki";
my $kfOut  = $out . ".kf";
my $lmxOut = $out . ".lmx";
my $zaOut  = $out . ".za";
print $crdOut,"\n";


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

my $NECP;
my $ECPatom;
my $ECPl;
my @ECPZa;
my @ECPList;
my @ECPKi;
my @ECPKf;
my @ECPLMax;
my @ECPA;
my @ECPa;
my $ECPNa;
my $ECPNk;

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
  }elsif($preLine =~ m/^Atomic numbers/){
    $preLine =~ s/.*N= *//;
    $ECPl = $preLine;
    unless($ECPl % 6){$Diff = 1;}
    else{$Diff = 0;}
    my $NECPLine = int($ECPl / 6) - $Diff;
    my $ECPtmp = $line;
    for(my $i=0;$i<$NECPLine;$i++){
      $ECPtmp .= <$fchk>;
    }
    @ECPZa = split(" ", $ECPtmp)
  }elsif($preLine =~ m/ECP-KFirst/){
    $preLine =~ s/.*N= *//;
    $ECPl = $preLine;
    unless($ECPl % 6){$Diff = 1;}
    else{$Diff = 0;}
    my $NECPLine = int($ECPl / 6) - $Diff;
    my $ECPtmp = $line;
    for(my $i=0;$i<$NECPLine;$i++){
      $ECPtmp .= <$fchk>;
    }
    @ECPKi = split(" ", $ECPtmp)
  }elsif($preLine =~ m/ECP-KLast/){
    $preLine =~ s/.*N= *//;
    $ECPl = $preLine;
    unless($ECPl % 6){$Diff = 1;}
    else{$Diff = 0;}
    my $NECPLine = int($ECPl / 6) - $Diff;
    my $ECPtmp = $line;
    for(my $i=0;$i<$NECPLine;$i++){
      $ECPtmp .= <$fchk>;
    }
    
    @ECPKf = split(" ", $ECPtmp)
  }elsif($preLine =~ m/ECP-LMax/){
    $preLine =~ s/.*N= *//;
    $ECPl = $preLine;
    unless($ECPl % 6){$Diff = 1;}
    else{$Diff = 0;}
    my $NECPLine = int($ECPl / 6) - $Diff;
    my $ECPtmp = $line;
    for(my $i=0;$i<$NECPLine;$i++){
      $ECPtmp .= <$fchk>;
    }
    @ECPLMax = split(" ", $ECPtmp)
  }elsif($preLine =~ m/ECP-LPSkip/){
    $preLine =~ s/.*N= *//;
    $NECP = $preLine;
    unless($NECP % 6){$Diff = 1;}
    else{$Diff = 0;}
    my $NECPLine = int($NECP / 6) - $Diff;
    $ECPatom = $line;
    for(my $i=0;$i<$NECPLine;$i++){
      $ECPatom .= <$fchk>;
    }
    @ECPList = split(" ", $ECPatom);
  }elsif($preLine =~ m/ECP-CLP1/){
    $preLine =~ s/.*N= *//;
    $ECPl = $preLine;
    unless($ECPl % 5){$Diff = 1;}
    else{$Diff = 0;}
    my $NECPLine = int($ECPl / 5) - $Diff;
    my $ECPtmp = $line;
    for(my $i=0;$i<$NECPLine;$i++){
      $ECPtmp .= <$fchk>;
    }
    @ECPA = split(" ", $ECPtmp)
  }elsif($preLine =~ m/ECP-ZLP/){
    $preLine =~ s/.*N= *//;
    $ECPl = $preLine;
    unless($ECPl % 5){$Diff = 1;}
    else{$Diff = 0;}
    my $NECPLine = int($ECPl / 5) - $Diff;
    my $ECPtmp = $line;
    for(my $i=0;$i<$NECPLine;$i++){
      $ECPtmp .= <$fchk>;
    }
    @ECPa = split(" ", $ECPtmp);
    $ECPNa = $NECPLine;
  }
  $preLine = $line;
}

foreach(@ECPZa){
  print $_,"\n";
}

# sum up number of AOs
my $sum = 0;
my $j=1;
my @AOMAP;
my @lmMAP;
for(my $i=0;$i<scalar(@NSHELL);$i++){
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
    case -2 {
      #print "$i\n",$sum+1,"-",$sum+6,"\n";
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
    case -3 {
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

my $c=0;


################################################################
#                                                              # 
#  File discriptions                                           #
#   1. CMNBasis.pltmp:                                         #
#       Charge,                                                #
#       Multiplicity,                                          #
#       Number of total e,                                     #
#       number ofalpha e,                                      #
#       number of beta e,                                      #
#       number of Basis functions,                             #
#       number of premitive gaussian functions,                #
#       number of contracted shells,                           #
#       number of atoms,                                       #
#      are printed in the above order                          #
#   2. TYPE.pltmp:                                             #
#       list of shell types, s=0, p=1, (cartisian)d= 2,        #
#                                      (spherical)d=-2         #
#   3. NSHELL.pltmp: Number of premitive gaussians per shell   #
#   4. AOMAP.pltmp: map from basis to contracted shell (AOs)   #
#   5. AGMAP.pltmp: map from shell to premitive gaussians      #
#   6. MAP.pltmp: map from contracted shell to atom            #
#                                                              #
#  ECP related files                                           #
#   7. ECPKi.pltmp: starting k index for each atom             #
#   8. ECPKf.pltmp: ending k index for each atom               #
#   9. ECPList.pltmp: list of activated ECP atom               #
#      (1: all-electron) (0: pp)                               #
#                                                              #
################################################################

open(my $fcoord, "> ./$crdOut") || die "open file\n";
for(my $i=0;$i<$NATOM;$i++){
  printf $fcoord "%d % 14.8E % 14.8E % 14.8E\n", 
  $ACHARGE[$i], $COORD[3*$i], $COORD[3*$i+1], $COORD[3*$i+2];
}
open(my $fecpki, "> ./$kiOut") || die "open file\n";
open(my $fecpkf, "> ./$kfOut") || die "open file\n";
open(my $fecplmax, "> ./$lmxOut") || die "open file\n";
open(my $fecplist, "> ./$pplOut") || die "open file\n";
open(my $fecpA, "> ./$alkOut") || die "open file\n";
open(my $fecpa, "> ./$expOut") || die "open file\n";
open(my $fecpza,"> ./$zaOut") || die "open file\n";
for(my $i=0;$i<scalar(@ECPKi);$i++){
  print $fecpki "$ECPKi[$i]\n";
  print $fecpkf "$ECPKf[$i]\n";
}
for(my $i=0;$i<scalar(@ECPLMax);$i++){
  print $fecplmax "$ECPLMax[$i]\n";
  print $fecplist "$ECPList[$i]\n";
  print $fecpza "$ECPZa[$i]\n";
}
for(my $i=0;$i<scalar(@ECPa);$i++){
  print $fecpa "$ECPa[$i]\n";
  print $fecpA "$ECPA[$i]\n";
}


#open(my $fcmn, "> ./CMNBasis.pltmp") || die "open file\n";
#open(my $facharge, "> ./ACHARGE.pltmp") || die "open file\n";
#open(my $ftype, "> ./TYPE.pltmp") || die "open file\n";
#open(my $fnshell, "> ./NSHELL.pltmp") || die "open file\n";
#open(my $fmap, "> ./MAP.pltmp") || die "open file\n";
#open(my $fcoord, "> ./COORD.pltmp") || die "open file\n";
#open(my $fpexp, "> ./PEXP.pltmp") || die "open file\n";
#open(my $fctrc, "> ./CTRC.pltmp") || die "open file\n";
#open(my $fewa, "> ./EW_a.pltmp") || die "open file\n";
#open(my $feva, "> ./EV_a.pltmp") || die "open file\n";
#open(my $faomap, "> ./AOMAP.pltmp") || die "open file\n";
#open(my $flmmap, "> ./lmMAP.pltmp") || die "open file\n";
#open(my $fagmap, "> ./AGMAP.pltmp") || die "open file\n";
#open(my $fecpki, "> ./ECPKi.pltmp") || die "open file\n";
#open(my $fecpkf, "> ./ECPKf.pltmp") || die "open file\n";
#open(my $fecplmax, "> ./ECPLMax.pltmp") || die "open file\n";
#open(my $fecplist, "> ./ECPList.pltmp") || die "open file\n";
#open(my $fecpA, "> ./ECPA.pltmp") || die "open file\n";
#open(my $fecpa, "> ./ECPa.pltmp") || die "open file\n";
#
#for(my $i=0;$i<scalar(@ECPKi);$i++){
#  print $fecpki "$ECPKi[$i]\n";
#  print $fecpkf "$ECPKf[$i]\n";
#}
#for(my $i=0;$i<scalar(@ECPLMax);$i++){
#  print $fecplmax "$ECPLMax[$i]\n";
#  print $fecplist "$ECPList[$i]\n";
#}
#for(my $i=0;$i<scalar(@ECPa);$i++){
#  print $fecpa "$ECPa[$i]\n";
#  print $fecpA "$ECPA[$i]\n";
#}
#
#
#foreach(@CMNBasis){
#  my $cmn = $_;
#  $cmn =~ s/\n//;
#  print $fcmn "$cmn\n";
#}
#print $fcmn scalar(@PEXP),"\n";
#print $fcmn scalar(@MAP),"\n";
#print $fcmn scalar(@ACHARGE),"\n";
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
#}
#for(my $i=0;$i<scalar(@AOMAP);$i++){
#  print $faomap "$AOMAP[$i]\n";
#  print $flmmap "$lmMAP[$i]\n";
#}
#for(my $i=0;$i<scalar(@EV_a);$i++){
#  print $feva "$EV_a[$i]\n";
#}
#for(my $i=0;$i<$NATOM;$i++){
#  print $fcoord "$COORD[3*$i] ";
#  print $fcoord "$COORD[3*$i+1] ";
#  print $fcoord "$COORD[3*$i+2]\n";
#  print $facharge "$ACHARGE[$i]\n";
#}
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
