#!/usr/bin/perl -w

use strict;
use warnings;

my $infile = join('',$ARGV[0]);
my $add = join('',$ARGV[1]);
my @nr;
my @hartree;
my @evolt;
my @signal;
my @numbers;
my $string;
my @superstring;
my @sorted_numbers;

open my $A, $infile or die "Could not open $infile:$!\n";
while (my $line = <$A>) {
    $string = $line;
    my $number = parseLine($line);
    push @superstring, $string;
    push @numbers, $number;
}
close($A);
@sorted_numbers = sort { $a <=> $b } @numbers;

foreach my $x (0..$#sorted_numbers){
    foreach my $y (0..$#superstring){
	if(index($superstring[$y],$sorted_numbers[$x]) != -1){
	    my $numm = parseNumm($superstring[$y]);
	    $numm += $add;
	    print "$numm: $superstring[$y]";
	}else{
	    next;
	}
    }
}


sub parseLine {
    my $str = shift;
    $str =~ /([\d]{1,5})[ \t]+([\d\.Ee\+\-]+)[ \t]+([\d\.]+)[ \t]+([\d\.Ee\+\-]+)/;
    my $stateno = $1;
    my $hartree = $2;
    my $evolt = $3;
    my $signal = $4;
    return $evolt;
}

sub parseNumm {
    my $str = shift;
    $str =~ /([\d]{1,5})[ \t]+([\d\.Ee\+\-]+)[ \t]+([\d\.]+)[ \t]+([\d\.Ee\+\-]+)/;
    my $stateno = $1;
    my $hartree = $2;
    my $evolt = $3;
    my $signal = $4;
    return $stateno;
}
