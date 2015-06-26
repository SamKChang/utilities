#!/usr/bin/perl

open(OUT,"<",$ARGV[0]) or die $!;
my $FLAG=0;
my $Et;


while(my $line=<OUT>){
  if($line =~ m/FINAL RESULTS/){
    $FLAG=1;
    while($line !~ m/======/){
      #print $line;
      if($line =~ m/TOTAL ENERGY/){
        $Et = $line;
        $Et =~ s/.*= {1,}//g;
        $Et =~ s/ A\.U\.//g;
        chop $Et;
        $FLAG += 1;
      }
      $line=<OUT>;
    }
  }
}

if($FLAG!=2){
  print "$ARGV[0]: file not complete...\n";
  die;
}


#print "K:$K\n";
#print "A:$A\n";
#print "S:$S\n";
#print "L:$L\n";
#print "N:$N\n";
#print "X:$X\n";
print "$Et\n";
