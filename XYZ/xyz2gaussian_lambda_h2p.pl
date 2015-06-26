#!/usr/bin/perl
use Cwd;
use Switch;
use List::Util qw(min max);
use POSIX;
use strict;
use warnings;

########################################################
# written for GEOALCHEMY project for X-H bond cpmdscan #
########################################################
# centering heavy atom at the YZ plan center
# YZ dimension and X length are specified as input
# xyz2cpmd_YZcenter.pl foo.xyz X YZ 

my $lambda = $ARGV[2];
my $nl = sprintf("%02d",$lambda*10);

my $pwd = cwd();

##################################
# read in number of atom $NA and #
# second line information $info2 #
##################################
open XYZ, "<", "$ARGV[0]" or die $!;
my $stem = $ARGV[0];
$stem =~ s/\..*//;
my $NA = <XYZ>;
$NA =~ s/\n//;
my $info2 = <XYZ>;
my $Nve = 0;
my $multi = 1;
my $charge = 0;
my @ATOMS;
my $itr;

##################################
# read atom type and coordinates #
##################################
for(my $itr=0;$itr<$NA;$itr++){
  # read and string into componets
  my @input = split(' ',<XYZ>);
  my @ATOM;
  $ATOM[0] = atom_a2n($input[0]);
  for(my $i=1;$i<4;$i++){$ATOM[$i] = $input[$i];}
  # assign multiD data array
  push (@{$ATOMS[$itr]}, @ATOM);
  $Nve += atom_ve($input[0]);
}
$multi += $Nve%2;
if($multi != 1){ 
  $charge = -1; 
  $multi = 1;
}
my $AStr = atom_a2n($ATOMS[1][0]);

###################
# second XYZ file #
###################
open XYZ2, "<", "$ARGV[1]" or die $!;
my $stem2 = $ARGV[1];
$stem2 =~ s/.*\/(.*)\..*/$1/;
$stem2 =~ s/\..*//;
my $NA2 = <XYZ2>;
my $info22 = <XYZ2>;
my $Nve2 = 0;
my $multi2 = 1;
my $charge2 = 0;
my @ATOMS2;

for(my $itr=0;$itr<$NA2;$itr++){
  # read and string into componets
  my @input = split(' ',<XYZ2>);
  my @ATOM;
  $ATOM[0] = atom_a2n($input[0]);
  for(my $i=1;$i<4;$i++){$ATOM[$i] = $input[$i];}
  # assign multiD data array
  push (@{$ATOMS2[$itr]}, @ATOM);
  $Nve2 += atom_ve($input[0]);
}
$multi2 += $Nve2%2;

# set multiplicity to 1 by adding electrons
if($multi2 != 1){
  $multi2 = 1;
  $charge2 = -1;
}
my $AStr2 = atom_a2n($ATOMS2[1][0]);

###########################################
# read what every info at the end of file #
###########################################
my @infoEND;
for($itr=0;<XYZ>;$itr++){
  #my @info = split(' ',$_);
  push (@infoEND, $_);
}
my $Ni=$itr;




#############################
# construct alchemical path #
#############################
my (@ATOMS_FIX, @ATOMS_A2B, @ATOMS_A2V, @ATOMS_V2B, @MAP_A2B);
@ATOMS_A2V = @ATOMS;
@ATOMS_V2B = @ATOMS2;
my ($n_fix, $n_a2b, $n_a2v, $n_v2b) = (0,0,0,0);
for(my $a1=0;$a1<$NA;$a1++){
  for(my $a2=0;$a2<$NA2;$a2++){
    my $R12=0;
    # check distance
    for(my $i=1;$i<4;$i++){
      $R12 += ($ATOMS[$a1][$i] - $ATOMS2[$a2][$i])**2;
    }
    $R12 = $R12**(1/2);
    # for close distance, check atom type
    if($R12 < 10**(-5)){
      my @data;
      $data[0] = $ATOMS[$a1][0];
      for(my $i=1;$i<4;$i++){
        $data[$i] = ($ATOMS[$a1][$i] + $ATOMS2[$a2][$i])/2.0
      }
      if($ATOMS[$a1][0] == $ATOMS2[$a2][0]){
        push(@{$ATOMS_FIX[$n_fix++]},@data);
      }else{
        push(@{$ATOMS_A2B[$n_a2b]},@data);
        $MAP_A2B[$n_a2b++] = $ATOMS2[$a2][0];
      }
      splice @ATOMS_A2V, $a1-($n_a2v++), 1;
      splice @ATOMS_V2B, $a2-($n_v2b++), 1;
    }
  }
}
$n_a2v = $NA - $n_a2v;
$n_v2b = $NA2 - $n_v2b;
## for debug purpose
## print coordinates for alchemical path
#print "ATOMS_FIX: $n_fix\n";
#for (my $a=0;$a<$n_fix;$a++){
#  for(my $i=0;$i<4;$i++){
#    print $ATOMS_FIX[$a][$i]," ";
#  }
#  print "\n";
#}
#print "\n";
#print "ATOMS_A2B: $n_a2b\n";
#for (my $a=0;$a<$n_a2b;$a++){
#  for(my $i=0;$i<4;$i++){
#    print $ATOMS_A2B[$a][$i]," ";
#  }
#  print "\n";
#}
#print "\n";
#print "ATOMS_A2V: $n_a2v\n";
#for (my $a=0;$a<$n_a2v;$a++){
#  for(my $i=0;$i<4;$i++){
#    print $ATOMS_A2V[$a][$i]," ";
#  }
#  print "\n";
#}
#print "\n";
#print "ATOMS_V2B: $n_v2b\n";
#for (my $a=0;$a<$n_v2b;$a++){
#  for(my $i=0;$i<4;$i++){
#    print $ATOMS_V2B[$a][$i]," ";
#  }
#  print "\n";
#}
#print "\n";


####################
## data processing #
####################
#my ($NTYPE_FIX,$NTYPE_A2B,$NTYPE_A2V,$NTYPE_V2B);
## fixed atoms
## numerical sort first column of @atoms
#@ATOMS_FIX = sort {$a->[0] <=> $b->[0]} @ATOMS_FIX;
#
## atom type map
#my @aMAP_FIX;
#$aMAP_FIX[0]=0;
#if($n_fix>0){
#  $NTYPE_FIX = 1;
#}else{
#  $NTYPE_FIX = 0;
#}
#for(my $a=1;$a<$n_fix;$a++){
#  if($ATOMS_FIX[$a][0] != $ATOMS_FIX[$a-1][0]){
#		$aMAP_FIX[$NTYPE_FIX++] = $a;
#  }
#}
#push(@aMAP_FIX,$n_fix);
#
#for(my $i=1;$i<4;$i++){
#  for(my $a=0;$a<$n_fix;$a++){
#    $ATOMS_FIX[$a][$i] += $CENTER;
#  }
#}
#
## mutating atoms
## numerical sort first column of @atoms
#@ATOMS_A2B = sort {$a->[0] <=> $b->[0]} @ATOMS_A2B;
#
## atom type map
#my @aMAP_A2B;
#$aMAP_A2B[0]=0;
#if($n_a2b>0){
#  $NTYPE_A2B = 1;
#}else{
#  $NTYPE_A2B = 0;
#}
#for(my $a=1;$a<$n_a2b;$a++){
#  if($ATOMS_A2B[$a][0] != $ATOMS_A2B[$a-1][0]){
#		$aMAP_A2B[$NTYPE_A2B++] = $a;
#  }
#}
#push(@aMAP_A2B,$n_a2b);
#
#for(my $i=1;$i<4;$i++){
#  for(my $a=0;$a<$n_a2b;$a++){
#    $ATOMS_A2B[$a][$i] += $CENTER;
#  }
#}
#
## vanishing atoms
## numerical sort first column of @atoms
#@ATOMS_A2V = sort {$a->[0] <=> $b->[0]} @ATOMS_A2V;
#
## atom type map
#my @aMAP_A2V;
#$aMAP_A2V[0]=0;
#if($n_a2v>0){
#  $NTYPE_A2V = 1;
#}else{
#  $NTYPE_A2V = 0;
#}
#for(my $a=1;$a<$n_a2v;$a++){
#  if($ATOMS_A2V[$a][0] != $ATOMS_A2V[$a-1][0]){
#		$aMAP_A2V[$NTYPE_A2V++] = $a;
#  }
#}
#push(@aMAP_A2V,$n_a2v);
#
#for(my $i=1;$i<4;$i++){
#  for(my $a=0;$a<$n_a2v;$a++){
#    $ATOMS_A2V[$a][$i] += $CENTER;
#  }
#}
#
## growing atoms
## numerical sort first column of @atoms
#@ATOMS_V2B = sort {$a->[0] <=> $b->[0]} @ATOMS_V2B;
#
## atom type map
#my @aMAP_V2B;
#$aMAP_V2B[0]=0;
#if($n_v2b>0){
#  $NTYPE_V2B = 1;
#}else{
#  $NTYPE_V2B = 0;
#}
#for(my $a=1;$a<$n_v2b;$a++){
#  if($ATOMS_V2B[$a][0] != $ATOMS_V2B[$a-1][0]){
#		$aMAP_V2B[$NTYPE_V2B++] = $a;
#  }
#}
#push(@aMAP_V2B,$n_v2b);
#
#for(my $i=1;$i<4;$i++){
#  for(my $a=0;$a<$n_v2b;$a++){
#    $ATOMS_V2B[$a][$i] += $CENTER;
#  }
#}
#
#
## for debug purpose
## print coordinates for alchemical path
#print "ATOMS_FIX: $n_fix\n";
#for (my $a=0;$a<$n_fix;$a++){
#  for(my $i=0;$i<4;$i++){
#    print $ATOMS_FIX[$a][$i]," ";
#  }
#  print "\n";
#}
#print "\n";
#print "ATOMS_A2B: $n_a2b\n";
#for (my $a=0;$a<$n_a2b;$a++){
#  for(my $i=0;$i<4;$i++){
#    print $ATOMS_A2B[$a][$i]," ";
#  }
#  print "\n";
#}
#print "\n";
#print "ATOMS_A2V: $n_a2v\n";
#for (my $a=0;$a<$n_a2v;$a++){
#  for(my $i=0;$i<4;$i++){
#    print $ATOMS_A2V[$a][$i]," ";
#  }
#  print "\n";
#}
#print "\n";
#print "ATOMS_V2B: $n_v2b\n";
#for (my $a=0;$a<$n_v2b;$a++){
#  for(my $i=0;$i<4;$i++){
#    print $ATOMS_V2B[$a][$i]," ";
#  }
#  print "\n";
#}
#print "\n";
#print "$NTYPE_A2B\n";













## shift to the corner of first quadrant
#for(my $i=1;$i<4;$i++){
#  my @vec = map { $_->[$i]} @ATOMS;
#  my $min = min @vec;
#  for(my $a=0;$a<$NA;$a++){
#    $ATOMS[$a][$i] -= $min;
#  }
#}
#
## shift molecule to MARGEIN from quadrant corner
#for(my $i=1;$i<4;$i++){
#  my @vec = map { $_->[$i]} @ATOMS;
#  my $max = max @vec;
#  $BOXD[$i-1] = $max + $PERIODIC_DISTANCE;
#}
#for(my $i=0;$i<3;$i++){
#  for(my $a=0;$a<$NA;$a++){
#    $ATOMS[$a][$i+1] += $margin;
#  }
#}

####################
# Gaussian09 NOTES #
####################
# command for Hamiltonian components print out
# IOp(3/5=17) request MWB default basis
# IOp(3/6=6) request MWB default ECP
# IOp(3/33=4) request print out 1-particle Matrix components
# IOp(2/12=3) supress "atoms too near" check
# guess(indo) request semiemperical initial guess for density matrix

#######################
# Gaussian09 commands #
#######################
# hf 6d 10f nosymm Scf(maxcycle=1,vshift,verytight) int(grid=ultrafine) IOp(3/5=17) IOp(3/6=6) IOp(3/33=4) IOp(2/12=3) guess(indo)
# ropbe1pbe 6d 10f nosymm Scf(maxcycle=1000,vshift,verytight) int(grid=ultrafine) IOp(3/5=17) IOp(3/6=6) IOp(2/12=3)
# ropbe1pbe 6d 10f nosymm Scf(maxcycle=1000,vshift,verytight) int(grid=ultrafine) IOp(3/5=17) IOp(3/6=6) IOp(3/33=4) guess(indo)
#%chk=$pwd/com/$stem
# ropbe1pbe 6d 10f nosymm Scf(maxcycle=1000,vshift,verytight) int(grid=superfine) IOp(3/5=17) IOp(3/6=6) IOp(3/33=4) guess(indo)
# pbepbe/def2-TZVP 6d 10f nosymm Scf(maxcycle=1000,verytight) int(grid=ultrafine)
#%chk=$pwd/$stem

###########################
# print Gaussian09 header #
###########################
print <<EOF;
%nproc=1
# pbepbe gen 6d 10f nosymm Scf(maxcycle=1000,verytight) int(grid=ultrafine) charge IOp(2/12=3)

$stem

$charge   $multi
EOF

###################################
# print atom type and coordinates #
###################################
for(my $a=0;$a<$n_a2v;$a++){
  printf("%2s   ", atom_a2n($ATOMS_A2V[$a][0]));
  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS_A2V[$a][$i]);}
  print "\n";
}
for(my $a=0;$a<$n_a2b;$a++){
  printf("%2s   ", atom_a2n($ATOMS_A2B[$a][0]));
  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS_A2B[$a][$i]);}
  print "\n";
}
for(my $a=0;$a<$n_v2b;$a++){
  printf("%2s-Bq", atom_a2n($ATOMS_V2B[$a][0]));
  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS_V2B[$a][$i]);}
  print "\n";
}
for(my $a=0;$a<$n_fix;$a++){
  printf("%2s   ", atom_a2n($ATOMS_FIX[$a][0]));
  for(my $i=1;$i<4;$i++){printf(" % 12.8f",$ATOMS_FIX[$a][$i]);}
  print "\n";
}

## print ghost atoms
#for(my $d=0.5;$d<=3.1;$d+=0.1){
#  if(abs($d-$ATOMS[0][1])>0.01 and abs($d-$ATOMS2[0][1])>0.01){
#    printf("%-2s", "H-Bq");
#    printf(" % 12.8f",$d);
#    printf(" % 12.8f",0);
#    printf(" % 12.8f",0);
#    print "\n";
#  }
#}

# print charges
print "\n";
for(my $a=0;$a<$n_a2v;$a++){
  for(my $i=1;$i<4;$i++){printf("% 12.8f ",$ATOMS_A2V[$a][$i]);}
  printf("% 4.2f\n",-$lambda)
}
for(my $a=0;$a<$n_a2b;$a++){
  for(my $i=1;$i<4;$i++){printf("% 12.8f ",$ATOMS_A2B[$a][$i]);}
  printf("% 4.2f\n",-$lambda*($ATOMS_A2B[$a][0]-$MAP_A2B[$a]));
}
for(my $a=0;$a<$n_v2b;$a++){
  for(my $i=1;$i<4;$i++){printf("% 12.8f ",$ATOMS_V2B[$a][$i]);}
  printf("% 4.2f\n",$lambda)
}



# print basis set
print <<EOF;

H     0 
S   1   1.00
     34.0613410              1.0000000
S   1   1.00
      5.1235746              1.0000000
S   1   1.00
      1.1646626              1.0000000
S   1   1.00
      0.32723041             1.0000000        
S   1   1.00
      0.10307241             1.0000000        
P   1   1.00
      0.8000000              1.0000000        
****


EOF


#####################################
## print atom type and coordinates #
####################################
## vanishing atoms
#for(my $t=0;$t<$NTYPE_A2V;$t++){
#  my $NAME_A = atom_a2n($ATOMS_A2V[$aMAP_A2V[$t]][0]);
#  my $NT = $aMAP_A2V[$t+1] - $aMAP_A2V[$t];
#  print "*$NAME_A","2V","_$ppFUNCTION","_",$nl,".psp ","FRAC\n";
#  print " LMAX=F\n";
#  print "  $NT\n";
#  for(my $a=$aMAP_A2V[$t];$a<$aMAP_A2V[$t+1];$a++){
#    print "  ";
#    for(my $i=1;$i<4;$i++){
#      printf("% 12.8f ",$ATOMS_A2V[$a][$i]);
#    }
#    print "\n";
#  }
#  print "\n";
#}
#
## mutating atoms
#for(my $t=0;$t<$NTYPE_A2B;$t++){
#  my $NAME_A = atom_a2n($ATOMS_A2B[$aMAP_A2B[$t]][0]);
#  my $NAME_B = atom_a2n($MAP_A2B[$t]);
#  my $NT = $aMAP_A2B[$t+1] - $aMAP_A2B[$t];
#  print "*$NAME_A","2","$NAME_B","_$ppFUNCTION","_",$nl,".psp ","FRAC\n";
#  print " LMAX=F\n";
#  print "  $NT\n";
#  for(my $a=$aMAP_A2B[$t];$a<$aMAP_A2B[$t+1];$a++){
#    print "  ";
#    for(my $i=1;$i<4;$i++){
#      printf("% 12.8f ",$ATOMS_A2B[$a][$i]);
#    }
#    print "\n";
#  }
#  print "\n";
#}
#
## growing atoms
#for(my $t=0;$t<$NTYPE_V2B;$t++){
#  my $NAME_A = atom_a2n($ATOMS_V2B[$aMAP_V2B[$t]][0]);
#  my $NT = $aMAP_V2B[$t+1] - $aMAP_V2B[$t];
#  print "*V2","$NAME_A","_$ppFUNCTION","_",$nl,".psp ","FRAC\n";
#  print " LMAX=F\n";
#  print "  $NT\n";
#  for(my $a=$aMAP_V2B[$t];$a<$aMAP_V2B[$t+1];$a++){
#    print "  ";
#    for(my $i=1;$i<4;$i++){
#      printf("% 12.8f ",$ATOMS_V2B[$a][$i]);
#    }
#    print "\n";
#  }
#  print "\n";
#}
#
## fixed atoms
#for(my $t=0;$t<$NTYPE_FIX;$t++){
#  my $NAME_A = atom_a2n($ATOMS_FIX[$aMAP_FIX[$t]][0]);
#  my $NT = $aMAP_FIX[$t+1] - $aMAP_FIX[$t];
#  print "*$NAME_A","_q",atom_ve("$NAME_A"),"_$ppFUNCTION",".psp ","FRAC\n";
#  print " LMAX=F\n";
#  print "  $NT\n";
#  for(my $a=$aMAP_FIX[$t];$a<$aMAP_FIX[$t+1];$a++){
#    print "  ";
#    for(my $i=1;$i<4;$i++){
#      printf("% 12.8f ",$ATOMS_FIX[$a][$i]);
#    }
#    print "\n";
#  }
#  print "\n";
#}
#
#
#####################
## print constraint #
#####################
##print <<EOF;
##CONSTRAINTS
## FIX ELEMENT
##  6
##END CONSTRAINTS
##EOF
#
#
#print "&END\n";
close(XYZ);

my @test;


###################
# useful routines #
###################

# count valence electrons
sub atom_a2n{
  switch($_[0]){
    case "H" {return  1}
    case "B" {return  5}
    case "C" {return  6}
    case "N" {return  7}
    case "O" {return  8}
    case "F" {return  9}
    case "Si"{return 14}
    case "P" {return 15}
    case "S" {return 16}
    case "Cl"{return 17}
    case "Ge"{return 32}
    case "As"{return 33}
    case "Se"{return 34}
    case "Br"{return 35}
    case "Sn"{return 50}
    case "Sb"{return 51}
    case "Te"{return 52}
    case "I" {return 53}
    case "1" {return "H" }
    case "5" {return "B" }
    case "6" {return "C" }
    case "7" {return "N" }
    case "8" {return "O" }
    case "9" {return "F" }
    case "14" {return "Si"}
    case "15" {return "P" }
    case "16" {return "S" }
    case "17" {return "Cl"}
    case "32" {return "Ge"}
    case "33" {return "As"}
    case "34" {return "Se"}
    case "35" {return "Br"}
    case "50" {return "Sn"}
    case "51" {return "Sb"}
    case "52" {return "Te"}
    case "53" {return "I" }
    #else    {
    #  print "element: $_[0] not found\n";
    #  die;
    #}
  }
}

sub atom_ve{
  switch($_[0]){
    case "H" {return 1}
    case "B" {return 3}
    case "C" {return 4}
    case "Si"{return 4}
    case "Ge"{return 4}
    case "Sn"{return 4}
    case "N" {return 5}
    case "P" {return 5}
    case "As"{return 5}
    case "Sb"{return 5}
    case "O" {return 6}
    case "S" {return 6}
    case "Se"{return 6}
    case "Te"{return 6}
    case "F" {return 7}
    case "Cl"{return 7}
    case "Br"{return 7}
    case "I" {return 7}
    else     {
      print "element: $_[0] not found\n";
      die;
    }
  }
}
