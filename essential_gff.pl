#!/usr/bin/env perl 
#
# v2 just uses TSV essentiality table generated in R 

use strict; 
use warnings; 

my $sample = shift @ARGV;   ## Sample name in the essentiality table 
my $tsv = shift @ARGV;      ## STR.ann_ess.tsv 
my $gff = shift @ARGV;      ## STR.tradis.gff 


open TSV,"<",$tsv or die "$!"; 
open GFF,"<",$gff or die "$!"; 

my $ESS = {}; 
my $index;                  ## column of our sample 

while (<TSV>) { 
  chomp;
  my @t = split /\t/; 
  if (m/ins_count/) { 
    for (my $i=0; $i < scalar @t; $i++) { 
      $index = $i if ($t[$i] eq $sample.".ins_count");
    } 
    print STDERR "Index for sample $sample is $index.\n"; 
  } else {
    my $id = $t[0]; 
    $ESS->{$id}->{ic} = $t[$index]; 
    $ESS->{$id}->{ii} = $t[$index+1]; 
    $ESS->{$id}->{type} = $t[$index+2];
  }  
} 

while (<GFF>) {
  chomp; 
  $_ =~ s/;$//g;  
  m/\tID=(.*?);/; 
  my $id = $1; 
  printf "%s;ins_index=%.5f(%d);tradis_type=%s\n",$_,$ESS->{$id}->{ii},$ESS->{$id}->{ic},$ESS->{$id}->{type};
}

close TSV; 
close GFF;  
