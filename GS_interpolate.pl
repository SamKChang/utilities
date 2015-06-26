#!/usr/bin/perl

#######################
# read PP1 parameters #
#######################

if("<$ARGV[0]" =~ /VOID/){
   $ZV1 = 0;
   for ($i=0;$i<5;$i++){$Rlc1[$i] = 0};
   for ($i=0;$i<3;$i++){$Rnl1[$i] = 0};
   for($n=0;$n<3;$n++){
      for($i=1;$i<=3;$i++){
         for($j=1;$j<=3;$j++){
            $H1[$n][$i][$j] = 0;
         }
      }
   }
}else{
   open(PP1, "<$ARGV[0]") or die "Couldn't open PP file, $!";
   for ($i=0;$i<5;$i++){$Rlc1[$i]=0};
   for ($i=0;$i<3;$i++){$Rnl1[$i]=0};
   for($n=0;$n<3;$n++){
      for($i=1;$i<=3;$i++){
         for($j=1;$j<=3;$j++){
            $H1[$n][$i][$j] = 0;
         }
      }
   }
   while(<PP1>){
      if($_ =~ /^  Z  =/){
         @lineArray = split ' ', $_;
         $Z = $lineArray[2];
      }
   
      if($_ =~ /^  ZV =/){
         @lineArray = split ' ', $_;
         $ZV1= $lineArray[2];
      }
   
      # loop for R_loc identification
      if($_ =~ /r\(loc\)/){
         # convert input data string to array
         @lineArray = split ' ', <PP1>;
         # define unified data format
         for($i=0;$i<=$#lineArray;$i++){
            $Rlc1[$i] = $lineArray[$i];
         }
      }
      # loop for R_nl identification
      if($_ =~ /r\([0-9]\)/){
         @I = split ' ', $_;
         $I[0] =~ s/[^0-9]//g;
         @lineArray = split ' ', <PP1>;
         $Rnl1[$I[0]] = $lineArray[0];
         for($i=1;$i<=$#lineArray;$i++){
            $H1[$I[0]][1][$i] = $lineArray[$i];
         }
         $nextline = <PP1>;
         $N = 2;
         while($nextline !~ /[a-z]/){
            @lineArray = split ' ', $nextline;
            for($i=0;$i<=$#lineArray;$i++){
               $H1[$I[0]][$N][$N+$i] = $lineArray[$i];
            }
            $N++;
            $nextline = <PP1>;
         }
         seek(PP1, -length($nextline), 1)
      }
   }
   close (PP1);
}
#######################
# read PP2 parameters #
#######################
if("<$ARGV[1]" =~ /VOID/){
   $ZV2 = 0;
   for ($i=0;$i<5;$i++){$Rlc2[$i] = 0};
   for ($i=0;$i<3;$i++){$Rnl2[$i] = 0};
   for($n=0;$n<3;$n++){
      for($i=1;$i<=3;$i++){
         for($j=1;$j<=3;$j++){
            $H2[$n][$i][$j] = 0;
         }
      }
   }
}else{
   open(PP2, "<$ARGV[1]") or die "Couldn't open PP file, $!";
   for ($i=0;$i<5;$i++){$Rlc2[$i]=0};
   for ($i=0;$i<3;$i++){$Rnl2[$i]=0};
   for($n=0;$n<3;$n++){
      for($i=1;$i<=3;$i++){
         for($j=1;$j<=3;$j++){
            $H2[$n][$i][$j] = 0;
         }
      }
   }
   while(<PP2>){
      if("<$ARGV[0]" =~ /VOID/){
         if($_ =~ /^  Z  =/){
            @lineArray = split ' ', $_;
            $Z = $lineArray[2];
         }
      }

      if($_ =~ /^  ZV =/){
         @lineArray = split ' ', $_;
         $ZV2= $lineArray[2];
      }
      # loop for R_loc identification
      if($_ =~ /r\(loc\)/){
         # convert input data string to array
         @lineArray = split ' ', <PP2>;
         # define unified data format
         for($i=0;$i<=$#lineArray;$i++){
            $Rlc2[$i] = $lineArray[$i];
         }
      }
      if($_ =~ /r\([0-9]\)/){
         @I = split ' ', $_;
         $I[0] =~ s/[^0-9]//g;
         @lineArray = split ' ', <PP2>;
         $Rnl2[$I[0]] = $lineArray[0];
         for($i=1;$i<=$#lineArray;$i++){
            $H2[$I[0]][1][$i] = $lineArray[$i];
         }
         $nextline = <PP2>;
         $N = 2;
         while($nextline !~ /[a-z]/){
            @lineArray = split ' ', $nextline;
            for($i=0;$i<=$#lineArray;$i++){
               $H2[$I[0]][$N][$N+$i] = $lineArray[$i];
            }
            $N++;
            $nextline = <PP2>;
         }
         seek(PP2, -length($nextline), 1)
      }
   }
   close (PP2);
}

##########################
# interpolate parameters #
##########################

$l = $ARGV[2];
if($ARGV[0] =~ m/VOID/){
  $Rlc[0] = $Rlc2[0];
  $Rnl[0] = $Rnl2[0];
}elsif($ARGV[1] =~ m/VOID/){
  $Rlc[0] = $Rlc1[0];
  $Rnl[0] = $Rnl1[0];
}else{
  $Rlc[0] = (1-$l)*$Rlc1[0] + $l*$Rlc2[0];
  $Rnl[0] = (1-$l)*$Rnl1[0] + $l*$Rnl2[0];
}
for ($i=1;$i<5;$i++){$Rlc[$i] = (1-$l)*$Rlc1[$i] + $l*$Rlc2[$i]};
$ZV = (1-$l)*$ZV1 + $l*$ZV2;
for ($i=1;$i<3;$i++){$Rnl[$i] = (1-$l)*$Rnl1[$i] + $l*$Rnl2[$i]};
for($n=0;$n<3;$n++){
   for($i=1;$i<=3;$i++){
      for($j=1;$j<=3;$j++){
         $h1= $H1[$n][$i][$j];
         $h2= $H2[$n][$i][$j];
         $H[$n][$i][$j] = (1-$l)*$h1 + $l*$h2;
      }
   }
}

##############
# print INFO #
##############

# local part
print "&INFO 
***********************Alchemical PPs***************************

 linear interpolated PPs 
 from PP1:$ARGV[0] to PP2:$ARGV[1] at lambda=$l

 local parameters:\n";
printf "       %10s %10s %10s %10s %10s\n", "r(loc)", "C(1)", "C(2)", "C(3)", "C(4)";
printf "  PP1: % 10f % 10f % 10f % 10f % 10f\n",$Rlc1[0],$Rlc1[1],$Rlc1[2],$Rlc1[3],$Rlc1[4];
printf "  PP2: % 10f % 10f % 10f % 10f % 10f\n\n",$Rlc2[0],$Rlc2[1],$Rlc2[2],$Rlc2[3],$Rlc2[4];
printf "  out: % 10f % 10f % 10f % 10f % 10f<---\n\n",$Rlc[0],$Rlc[1],$Rlc[2],$Rlc[3],$Rlc[4];

# nonlocal part r(0)
printf " Nonlocal parameters:\n";
printf "       %10s %10s\n","r(0)","h(i,j)^0";
printf "  PP1: % 10f % 10f % 10f % 10f\n",$Rnl1[0],$H1[0][1][1],$H1[0][1][2],$H1[0][1][3];
printf "       %10s %10s % 10f % 10f\n","","",$H1[0][2][2],$H1[0][2][3];
printf "       %10s %10s %10s % 10f\n","","","",$H1[0][3][3];
printf "  PP2: % 10f % 10f % 10f % 10f\n",$Rnl2[0],$H2[0][1][1],$H2[0][1][2],$H2[0][1][3];
printf "       %10s %10s % 10f % 10f\n","","",$H2[0][2][2],$H2[0][2][3];
printf "       %10s %10s %10s % 10f\n\n","","","",$H2[0][3][3];
printf "  out: % 10f % 10f % 10f % 10f <---\n",$Rnl[0],$H[0][1][1],$H[0][1][2],$H[0][1][3];
printf "       %10s %10s % 10f % 10f <---\n","","",$H[0][2][2],$H[0][2][3];
printf "       %10s %10s %10s % 10f <---\n\n","","","",$H[0][3][3];

# nonlocal part r(1)
printf "       %10s %10s\n","r(1)","h(i,j)^1";
printf "  PP1: % 10f % 10f % 10f\n",$Rnl1[1],$H1[1][1][1],$H1[1][1][2];
printf "       %10s %10s % 10f\n","","",$H1[1][2][2];
printf "  PP2: % 10f % 10f % 10f\n",$Rnl2[1],$H2[1][1][1],$H2[1][1][2];
printf "       %10s %10s % 10f\n\n","","",$H2[1][2][2];
printf "  out: % 10f % 10f % 10f <---\n",$Rnl[1],$H[1][1][1],$H[1][1][2];
printf "       %10s %10s % 10f <---\n\n","","",$H[1][2][2];
printf "       %10s %10s\n","r(2)","h(i,j)^2";
printf "  PP1: %10f % 10f\n",$Rnl1[2],$H1[2][1][1];
printf "  PP2: %10f % 10f\n\n",$Rnl2[2],$H2[2][1][1];
printf "  out: %10f % 10f <---\n\n",$Rnl[2],$H[2][1][1];
print"***********************Alchemical PPs***************************";

####################
# print PPs format #
####################

print"
&END
&ATOM
 Z  = $Z
 ZV = $ZV
 XC = 1134   0.6666666667
 TYPE = NORMCONSERVING GOEDECKER
&END
&POTENTIAL
    GOEDECKER
3                                    LMAX
  $Rlc[0]                                     RC
  4 $Rlc[1] $Rlc[2] $Rlc[3] $Rlc[4] #C C1 C2 C3 C4
  $Rnl[0] 3 $H[0][1][1] $H[0][1][2] $H[0][1][3] $H[0][2][2] $H[0][2][3] $H[0][3][3] H(s) 11 12 13 22 23 33
  $Rnl[1] 2 $H[1][1][1] $H[1][1][2] $H[1][2][2] H(p) 11 12 22
  $Rnl[2] 1 $H[2][1][1] H(d) 11
&END\n";
