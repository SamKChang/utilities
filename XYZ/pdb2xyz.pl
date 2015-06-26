#!/usr/bin/perl

my @xyz;
my $itr = 0;
open(PDB,$ARGV[0]);

while(my $line=<PDB>){
  chop($line);
  $line =~ s/ {1,}/ /g;
  my @data = split(" ", $line);
  if(($data[0] =~ m/ATOM/)||($data[0] =~ m/HETATM/)){
    $data[scalar(@data)-1] =~ s/[0-9].*//;
    if(scalar(@data)==12){
      $xyz[$itr][0] = $data[11];
      $xyz[$itr][1] = $data[6];
      $xyz[$itr][2] = $data[7];
      $xyz[$itr][3] = $data[8];
    }elsif(scalar(@data==11)){
      $xyz[$itr][0] = $data[10];
      $xyz[$itr][1] = $data[5];
      $xyz[$itr][2] = $data[6];
      $xyz[$itr][3] = $data[7];
    }elsif(scalar(@data==10)){
      $xyz[$itr][0] = $data[9];
      $xyz[$itr][1] = $data[4];
      $xyz[$itr][2] = $data[5];
      $xyz[$itr][3] = $data[6];
    }else{
      print "DATA FORMAT NOTMATCHED!\n";
      exit;
    }
    $itr += 1;
  }
}

print $itr,"\n\n";
for(my $i=0;$i<$itr;$i++){
  printf("%-2s",$xyz[$i][0]);
  for(my $j=0;$j<3;$j++){
    printf("% 14.8f",$xyz[$i][$j+1]);
  }
  print "\n";
}
close(PDB);
