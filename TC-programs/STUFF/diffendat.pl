#!/usr/bin/perl -w

use strict;
use warnings;

my $infile = join('',$ARGV[0]);
my $lineno = 0;
my @string;
my $number;
my $oldnumber;
my $energy;
my $oldenergy;
my $diff;

open my $A, $infile or die "Could not open $infile:$!\n";
while(my $line = <$A>) {
    if($lineno > 0){
	$oldnumber = $number;
	$oldenergy = $energy;
    }
    @string = parseline($line);
    $number = $string[0];
    $energy = $string[1];
    $diff = $energy - $oldenergy;
    if($lineno == 0){
	printf("%.4d    %12.10f\n",$number,$energy);
    }else{
	printf("%.4d    %12.10f   %12.10f\n",$number,$energy,$diff);
    }
    $lineno++;
}
close($A);

sub parseline {
    my $str = shift;
    $str =~ /([\d]{1,4})[ \t]+([\d\.\+\-Ee]+)/;
    my @num;
    $num[0] = $1;
    $num[1] = $2;
    return @num;
}

