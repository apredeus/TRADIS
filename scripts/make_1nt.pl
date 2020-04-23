#!/usr/bin/env perl 

use strict; 
use warnings; 

while (<>) { 
  if (m/^@/) { 
    ## simply print if header
    print; 
    next; 
  }

  chomp;  
  my @t = split /\t/; 
  ## following processing is clumsy but avoids special character mess that happens in qual string 
 
  my $read = shift @t;
  my $flag = shift @t; 
  my $chr  = shift @t; 
  my $pos  = shift @t; 
  my $mapq = shift @t; 
  my $cigar = shift @t;
  my $rn   = shift @t; 
  my $pn   = shift @t; 
  my $tlen = shift @t;
  my $seq  = shift @t;
  my $qual = shift @t; 
  
  my $rest = join "\t",@t; 

  if ($flag&16) { 
    ## things get quite tricky when you are on (-) strand
    $flag = 16; 
    $pos = $pos + ref_length($cigar) - 1; 
    $seq = substr($seq,-1,1); 
    $qual = substr($qual,-1,1); 
  } else { 
    $flag = 0; 
    $seq = substr($seq,0,1);
    $qual = substr($qual,0,1);
  } 
  # well..  
  print "$read\t$flag\t$chr\t$pos\t$mapq\t1M\t*\t0\t0\t$seq\t$qual\t$rest\n";
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
