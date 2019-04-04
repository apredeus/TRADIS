#!/usr/bin/perl
use strictures 1;
LINE: while (<>) {
    my @F=split /\t/, $_;
    if ($F[0]=~/^@/) {print; next}
    my @reads;
    push @reads, [@F];
    process_group(@reads) if eof;
    while (<>) {
        my @F=split /\t/, $_;
        if (not $F[0] eq $reads[0][0]) {
            process_group(@reads);
            redo LINE;
        } else {
            push @reads, [@F];
            process_group(@reads) if eof;
        }
    }
}

sub process_group {
    my @reads=@_;
    my @print;
    for my $r (@reads) {
        my @F=@$r;
        if ($F[1]&64 and (($F[1]&16 and $F[5]=~/\d+M$/) or ($F[1]^16 and $F[5]=~/^\d+M/))) {
            push @print, 1;
        } elsif ($F[1]&128) {            
            $r->[1]+=256 if $F[1]&2048;
            push @print, 1;
        } else {
            push @print, 0;
        }
#        say STDERR join "\t", @F, $print[-1];
    }
    if (not grep {$_==0} @print) {
        for my $r (@reads) {
            next if $r->[9] eq '*';
            next if $r->[9] eq '';
            print join "\t", @$r
        }
    }
} 
