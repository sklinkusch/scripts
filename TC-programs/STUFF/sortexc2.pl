#!/usr/bin/perl -w

use strict;
use warnings;

my $file = join('',$ARGV[0]);
my $no   = join('',$ARGV[1]);
my @string;
my $curr_state;
my $curr_frommo;
my $curr_tomo;
my $curr_val;
my %coeff;

open my $log, $file or die "Could not open $file: $!";
while (my $line = <$log>){
    if($line =~ /[\d]{1,4}[ \t]+[\d]{1,4}[ \t]+[\d]{1,4}[ \t]+[-]*[\d]\.[\d]{1,20}/){
	@string = parseline($line);
	$curr_state = $string[0];
	$curr_frommo = $string[1];
	$curr_tomo = $string[2];
	$curr_val = $string[3];
	if($curr_state == $no){
	    my $key = "$curr_frommo -> $curr_tomo";
	    $coeff{$key} = $curr_val;
	}elsif($curr_state > $no){
	    last;
	}
    }
}
close $log;

foreach my $key (sort { $coeff{$b} <=> $coeff{$a} } keys %coeff){
    printf "%12s: %10.8f\n", $key, $coeff{$key};
}

sub parseline {
    my $str = shift;
    $str =~ /([\d]{1,4})[ \t]+([\d]{1,4})[ \t]+([\d]{1,4})[ \t]+[-]*([\d\.]{3,22})/;
    my @numm;
    $numm[0] = $1;
    $numm[1] = $2;
    $numm[2] = $3;
    $numm[3] = $4;
    return @numm;
}

