#!/usr/bin/env perl 

use strict; 
use warnings; 
use Data::Dumper; 

my $param = shift @ARGV; 
my $gff = shift @ARGV; 

my $essen; 
my $ambig; 
my $genes = {}; 

my $csv = $param; 
my $tag = $param; 
$csv =~ s/.param.out//g; 
$tag =~ s/.tag.out.gz.*//g;

open PARAM,"<",$param or die "$!"; 
open GFF,"<",$gff or die "$!"; 
open CSV,"<",$csv or die "$!"; 

while (<PARAM>) { 
  chomp; 
  if (! m/essen/) {
    $essen = (split/\t/)[0]; 
    $ambig = (split/\t/)[1];
  } 
} 

print STDERR "Processing sample $tag; essentiality cutoff is $essen, ambiguous cutoff is $ambig\n"; 

<CSV>; 
while (<CSV>) { 
  chomp; 
  my @t = split /\t/;

  my $id = ($t[0] =~ m/([A-Z0-9]+_\d+)_?/) ? $1 : $t[0]; ## biotradis does weird shit with tRNA/rRNA locus tags 
  my $type = "non-essential"; 

  if ($t[7] < $essen) { 
    $type = "essential"; 
  } elsif ($t[7] < $ambig) { 
    $type = "ambiguous"; 
  }
 
  $genes->{$id}->{ii} = $t[7]; 
  $genes->{$id}->{type} = $type;
}

#print Dumper $genes; 

while (<GFF>) {
  chomp; 
  $_ =~ s/;$//g;  
  m/\tID=(.*?);/; 
  my $id = $1; 
#  printf STDERR "DEBUG: %s %s %s\n",$id,$genes->{$id}->{ii},$genes->{$id}->{type};
  printf "%s;ins_index=%.5f;tradis_type=%s\n",$_,$genes->{$id}->{ii},$genes->{$id}->{type}; 
}

close PARAM; 
close GFF; 
close CSV;  
