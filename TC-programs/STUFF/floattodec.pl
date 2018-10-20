#!/usr/bin/perl -w

use strict;
use warnings;

my $file = join('',$ARGV[0]);

open my $A, $file or die "Could not open $file: $!\n";
while(my $line = <$A>){
    my @string = split(/ /,$line);
    printf("%3u   %10.8f   %10.8f  % 12.10f\n",$string[0],$string[1],$string[2],$string[3]);
}

