#!/usr/bin/perl -w

use strict;
use warnings;

my $filename;
my @pops;
my $char = join('',$ARGV[0]);
my @enpop_x = (4,10,20,50,62,100,150,186,297);
#my $enpop_x = 297;

foreach my $x (0..$#enpop_x){
    $filename = sprintf("%spi-test-3pi%u.sts", $char, $enpop_x[$x]);
    open(my $A, "tail -1 $filename | ");
    while(my $line = <$A>){
	chomp $line;
	@pops = split(/[ ]+/,$line);
	my $fin = $#pops;
#	if($x != $#enpop_x){
	    printf("%3u   %22.20f   %22.20f   %22.20f\n", \
	    $enpop_x[$x], $pops[2], $pops[4], $pops[$fin]);
#	print "$enpop_x[$x]  $pops[0]  $pops[1]   $pops[2]    $pops[3]   $pops[4]   $pops[$fin-1]\n";
#    }else{
#	print "$enpop_x[$x]  $pops[0]  $pops[1]   $pops[2]    $pops[3]   $pops[4]   $pops[$fin]\n";
#    }
#	foreach my $y (0..$fin){
#	    printf("%3u   %22.20f\n", $y, $pops[$y]);
	    shift @pops;
#	}
	last;
    }
    close($A);
}

