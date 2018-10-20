#!/usr/bin/perl -w

use strict;
use warnings;

if($#ARGV == -1){
 print "Usage: ./average.pl <value1> <value2> ...\n";
 print "At least one value required!\n";
 exit;
}

my $num = 0;
my $sum = 0;
while($#ARGV > -1){
 $sum += $ARGV[0];
 $num += 1;
 shift(@ARGV);
}

my $res = $sum / $num;
print "Result: $res\n";

