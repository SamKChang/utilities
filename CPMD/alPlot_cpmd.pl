#!/usr/bin/perl

open(OUT,"<",$ARGV[0]) or die $!;
my $FLAG=0;
my ($K,$A,$S,$L,$N,$X);


while(my $line=<OUT>){
  if($line =~ m/FINAL RESULTS/){
    $FLAG=1;
    while($line !~ m/======/){
      #print $line;
      if($line =~ m/\(K\)/){
        $K = $line;
        $K =~ s/.*= {1,}//g;
        $K =~ s/ A\.U\.//g;
        chop $K;
        $FLAG += 1;
      }elsif($line =~ m/\(A\)/){
        $A = $line;
        $A =~ s/.*= {1,}//g;
        $A =~ s/ A\.U\.//g;
        chop $A;
        $FLAG += 1;
      }elsif($line =~ m/\(S\)/){
        $S = $line;
        $S =~ s/.*= {1,}//g;
        $S =~ s/ A\.U\.//g;
        chop $S;
        $FLAG += 1;
      }elsif($line =~ m/\(L\)/){
        $L = $line;
        $L =~ s/.*= {1,}//g;
        $L =~ s/ A\.U\.//g;
        chop $L;
        $FLAG += 1;
      }elsif($line =~ m/\(N\)/){
        $N = $line;
        $N =~ s/.*= {1,}//g;
        $N =~ s/ A\.U\.//g;
        chop $N;
        $FLAG += 1;
      }elsif($line =~ m/\(X\)/){
        $X = $line;
        $X =~ s/.*= {1,}//g;
        $X =~ s/ A\.U\.//g;
        chop $X;
        $FLAG += 1;
      }
      $line=<OUT>;
    }
  }
}

if($FLAG!=7){
  print "$ARGV[0]: file not complete...\n";
  die;
}


#print "K:$K\n";
#print "A:$A\n";
#print "S:$S\n";
#print "L:$L\n";
#print "N:$N\n";
#print "X:$X\n";
my $E = $K+$A-$S+$L+$N+$X;
print "$E\n";
