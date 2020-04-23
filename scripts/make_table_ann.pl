#!/usr/bin/env perl 
##
## simple script to convert GFF file to 
## number of lines in the output = number of unique IDs  

use strict; 
use warnings; 

my $gff = shift @ARGV; 
open GFF,"<",$gff or die "$!"; 

print "Locus_tag\tName\tChr\tBegin\tEnd\tStrand\tProduct\n";

while (<GFF>) { 
  chomp; 
  my @t = split /\t/; 
  if ($t[8] =~ m/ID=(.*?);/) {
    my $id = $1; 
    my $name = ($t[8] =~ m/Name=(.*?);/) ? $1 : "NONE"; 
    my $prod = ($t[8] =~ m/note=(.*?);/) ? $1 : "NONE";
    $prod =~ s/%2C/,/g;  
    print "$id\t$name\t$t[0]\t$t[3]\t$t[4]\t$t[6]\t$prod\n";
  } else { 
    next; 
  } 
}

close GFF;  
