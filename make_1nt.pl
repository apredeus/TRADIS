#!/usr/bin/env perl 

use strict; 
use warnings; 

while (<>) { 
  if (m/^@/) { 
    print; 
    next; 
  } 

  my @t = split /\t/; 
  s/\t$t[5]\t$t[6]\t$t[7]\t$t[8]\t/\t1M\t*\t0\t0\t/; 
  
  my $flag = $t[1]; 
  my $chr = $t[2]; 
  my $pos = $t[3];
  my $cigar = $t[5]; 
  my $seq = $t[9]; 
  my $qual = $t[10]; 

  if ($flag&16) { 
    ## things get quite tricky when you are on (-) strand
    $flag = 16; 
    $pos = $pos + ref_length($cigar) - 1; 
    $seq = substr($seq,-1,1); 
    $qual = substr($qual,-1,1); 
    s/\t$t[1]\t$t[2]\t$t[3]\t/\t$flag\t$chr\t$pos\t/; 
    s/\t$t[9]\t$t[10]\t/\t$seq\t$qual\t/; 
  } else { 
    $flag = 0; 
    $seq = substr($seq,0,1);
    $qual = substr($qual,0,1);
    s/\t$t[1]\t$t[2]\t$t[3]\t/\t$flag\t$chr\t$pos\t/; 
    s/\t$t[9]\t$t[10]\t/\t$seq\t$qual\t/;
  } 

  print; 
} 

sub ref_length {
  my $cigar = $_[0];
  my $len = 0;
  while ($cigar ne "") {
    $cigar =~ s/^(\d+)(\D+)//;
    $len += $1 if ($2 eq "M" || $2 eq "D" || $2 eq "N" || $2 eq "=" || $2 eq "X" || $2 eq "P");
  }
  return($len);
}
