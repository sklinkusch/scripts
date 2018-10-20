#!/usr/bin/perl -w

use strict;
use warnings;

# calculate percentage of margin

if($#ARGV != 1){
 print "Usage: ./margin.pl <max> <min>\n";
 exit;
}
my $max;
my $min;
if($ARGV[0] >= $ARGV[1]){
 $max = $ARGV[0];
 $min = $ARGV[1];
}else{
 $max = $ARGV[1];
 $min = $ARGV[0];
}

my $result = ($max - $min) / ($max + $min);
$result *= 100;

print "$result\n";
