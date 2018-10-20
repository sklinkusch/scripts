#!/usr/bin/perl -w

use strict;
use warnings;

my $file = join('',$ARGV[0]);

open my $A, $file or die "Could not open $file: $!";
while( my $line = <$A>) {
    $line =~ s/^[ ]+//;
    $line =~ s/\n//;
    my @string = split(/[ ]{2,7}/,$line);
    my $time = $string[0];
    my $timefs = $time/41.341373;
    shift @string;
    my @band;
    foreach my $x(0..$#string){
	if($x == 0){
	    $band[0] = $string[$x];             # ground state
	}elsif($x == 1){
	    $band[1] += $string[$x];            # first excited state
	}elsif($x == 2){
	    $band[2] += $string[$x];            # second excited state
	}elsif($x == 3){
	    $band[3] += $string[$x];            # third excited state
	}elsif($x == 4){
	    $band[4] += $string[$x];            # fourth excited state
	}elsif($x == 5){
	    $band[5] += $string[$x];            # fifth excited state
	}elsif($x == 6){
	    $band[6] += $string[$x];            # sixth excited state
	}elsif($x > 6 and $x < $#string){
	    $band[7] += $string[$x];            # higher bound states
	}elsif($x == $#string){
	    $band[8] += $string[$x];            # ion yield
	}
    }
    printf("%9.4f  %9.4f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n", $time,$timefs,$band[0],$band[1],$band[2],$band[3],$band[4],$band[5],$band[6],$band[7],$band[8]);
}
close($A);
