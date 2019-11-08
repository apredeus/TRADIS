#!/usr/bin/env perl
# this script removes all reads that do not map right after TRADIS tag
# it assumes that only R1 (-f 64) can have TRADIS tag on it  

use strict; 
use warnings; 

while (<>) {
  ## keep the header
  if (m/^@/) { 
    print;
    next;
  }
  ## assume only R1  
  my @t = split /\t/;
  print if ( ($t[1]&16 && $t[5] =~ m/\d+M$/) || ( !($t[1]&16) && $t[5] =~ m/^\d+M/)); 
} 
